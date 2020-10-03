function o = sample(problem, N, opts)
%Input:
% problem - a struct containing the constraints defining the polytope
%
% N - number of indepdent samples
% opts - sampling options
%
%Output:
% o - a structure containing the following properties:
%   samples - a dim x N vector containing all the samples
%   sampler - the subroutine to do one step
%   mixingTime - number of records per "independent" sample based on min(ess)
%   prepareTime - time to pre-process the input (including find analytic
%   center, remove redundant constraints, reduce dimension etc.)
%   sampleTime - total sampling time in seconds

t = tic;

%% Initialize parameters
o = struct;
o.iterPerRecord = 1;
o.checkESSIter = opts.checkESSIter;
o.samples = [];
o.stepSize = opts.initalStepSize;
o.momentum = 1 - min(1, o.stepSize/opts.effectiveStepSize);
o.iterSinceShrink = 0;
o.miniStepSinceShrink = 0;
o.rejectSinceShrink = 0;
o.done = false;

% set random seed
rng(opts.seed, 'simdTwister');
r = rng;
o.seed = r.Seed;

%% Presolve
P = Polytope(problem, opts);
ham = Hamiltonian(P, opts);
o.prepareTime = toc(t);
o.dim = ham.n - size(P.A,1);
o.ham = ham;
o.nSamples = 0;

%% Presolve
o.consecutiveBadStep = 0;
o.i = 1;
o.x = P.center;
o.v = o.ham.resample(o.x, zeros(size(o.x)), 0);
o.acceptedSteps = 0;

if opts.debugMode
    o.averageAccuracy = 0;
end

if opts.dynamicBound
    o.bound = max(abs(o.x), 1);
    ham.barrier.extraHessian = 1./(o.ham.n*o.bound.*o.bound);
end

if opts.weightedBarrier
    o.averageLSC = zeros(size(o.x));
end

ham.opts.checkPrecision = true;
while true
    % v step
    o.v = ham.resample(o.x, o.v, o.momentum);
    
    % x step
    o.x_last = o.x; o.v_last = o.v;
    o.H1 = ham.H(o.x, o.v);
    [o.x, o.v, o.ODEStep] = opts.odeMethod(o.x, o.v, o.stepSize, ham, opts);
    
    if ~ham.feasible(o.x, o.v)
        o.prob = 0;
    else
        o.H2 = ham.H(o.x, -o.v);
        o.prob = min(1, exp(o.H1-o.H2));
        
        if opts.dynamicBound
            o.bound = max(o.bound, abs(o.x));
        end
    end
    
    % rejection
    if rand >= o.prob
        opts.outputFunc('sample:reject', 'iter %i: rejected\n', o.i);
        o.x = o.x_last; o.v = -o.v_last; % flip the direction
        o.accepted = false;
    else
        o.accepted = true;
        o.acceptedSteps = o.acceptedSteps + 1;
    end
    ham.opts.checkPrecision = ~o.accepted; % check solver precision if not accepted
    
    % implement the dynamic bound
    if opts.dynamicBound && o.nSamples < opts.freezeMCMCAfterSamples
        idx = find((1./(o.bound.*o.bound) < o.ham.n*ham.barrier.extraHessian));
        if ~isempty(idx)
            opts.outputFunc('sample:bound_change', 'iter %i: changed the bound of %i coordinates.\n', o.i, length(idx));
            ham.barrier.extraHessian = 0.1./(o.ham.n*o.bound.*o.bound);
            o.v = o.ham.resample(o.x, zeros(size(o.x)), 0);
        end
    end
    
    if opts.weightedBarrier && o.nSamples < opts.freezeMCMCAfterSamples
        o.averageLSC = 0.9 * o.averageLSC + 0.1 * ham.last_lsc;
        idx = find(o.averageLSC > (0.5 / o.stepSize) * ham.barrier.weight);
        if ~isempty(idx)
            opts.outputFunc('sample:weight_change', 'iter %i: changed the weight of %i coordinates.\n', o.i, length(idx));
            ham.barrier.weight(idx) = min(1, ham.barrier.weight(idx) * 2);
            o.averageLSC = zeros(size(o.x));
        end
    end
    
    % if ODE solver takes too many step, check if higher accuracy needed
    if o.ODEStep > 2 * opts.targetODEStep
        ham.opts.checkPrecision = true;
    end
    
    if o.prob < 0.5 || o.ODEStep == opts.maxODEStep
        o.consecutiveBadStep = o.consecutiveBadStep + 1;
    else
        o.consecutiveBadStep = 0;
    end
    
    % Shrink the step during the warmup phase
    if opts.adaptiveStepSize && o.nSamples < opts.freezeMCMCAfterSamples
        shrink = 0;
        o.shiftedIter = o.iterSinceShrink + 20 / (1-o.momentum);
        o.rejectSinceShrink = o.rejectSinceShrink + 1-o.prob;
        
        targetProbability = (1-o.momentum)^(2/3) / 2;
        if (o.rejectSinceShrink > targetProbability  * o.shiftedIter)
            shrink = 'failure probability';
        end
        
        if (o.consecutiveBadStep > opts.consecutiveBadStepLimit)
            shrink = 'consecutive bad steps';
        end
        
        o.miniStepSinceShrink = o.miniStepSinceShrink + o.ODEStep;
        if (o.miniStepSinceShrink > opts.targetODEStep * o.shiftedIter)
            shrink = 'ODE accuracy';
        end
        
        if ischar(shrink)
            o.x = o.x_last; o.v = o.v_last;
            o.rejectSinceShrink = 0;
            o.miniStepSinceShrink = 0;
            o.iterSinceShrink = 0;
            o.stepSize = o.stepSize / opts.shrinkFactor;
            o.momentum = 1 - min(1, o.stepSize/opts.effectiveStepSize);
            opts.outputFunc('sample:step_adjust', 'Iter %i: Step shrinks due to %s. h=%f \n', o.i, shrink, o.stepSize);

            if o.stepSize < opts.minStepSize
                warning('iter %i: Algorithm fails to converge even with step size h=%f.\n', o.i, o.stepSize);
                save
                break;
            end
            o.consecutiveBadStep = 0;
        end
        
        o.iterSinceShrink = o.iterSinceShrink + 1;
    elseif o.consecutiveBadStep > opts.consecutiveBadStepLimit && o.nSamples >= opts.freezeMCMCAfterSamples
         % restart if too many bad steps after the warm up phase
        idx = randi(size(o.samples,2));
        o.x = o.samples(:,idx);
        o.v = o.ham.resample(o.x, zeros(size(o.x)), 0);
        o.iterSinceShrink = 0;
        o.miniStepSinceShrink = 0;
        o.rejectSinceShrink = 0;
        o.consecutiveBadStep = 0;
    end
    
    elapsed = toc(t);
    if (elapsed > opts.maxTime)
        warning('iter %i: opts.maxTime (%f sec) reached.\n', o.i, opts.maxTime);
        break;
    end
    
    if mod(o.i, o.iterPerRecord) == 0
        o.samples(:,end+1) = o.x;
    end
    
    if (o.i >= o.checkESSIter && o.acceptedSteps > opts.minAcceptedSteps)
        len = size(o.samples,2);
        o.ess = effective_sample_size(o.samples);
        o.nSamples = min(o.ess);
        o.mixingTime = o.i / o.nSamples;
        
        % thinning
        if opts.recordsPerIndependentSample < Inf
            h = max((len / o.nSamples) / opts.recordsPerIndependentSample, 1);
            idx = ceil(h*(1:floor(len / h))) + (len - ceil(h * floor(len / h)));
            o.samples = o.samples(:, idx);
            o.iterPerRecord = ceil(o.mixingTime / (opts.recordsPerIndependentSample * 2));
            if o.nSamples < N/2
                o.iterPerRecord = ceil(o.iterPerRecord / 2);
            end
        end
        
        if o.nSamples > N
            o.done = true;
            break;
        end
        
        % check ess after "freezeMCMCAfterSamples" samples and N/2 samples
        N_target = N;
        if o.nSamples < N/2
            N_target = min(N_target, N/2);
        end
        
        if o.nSamples < opts.freezeMCMCAfterSamples
            N_target = min(N_target, opts.freezeMCMCAfterSamples);
        end
        
        o.checkESSIter = o.i + (N_target - o.nSamples) * o.mixingTime;
    end
    
    if o.i > opts.maxStep
        warning('iter %i: opts.maxStep reached.\n', o.i, opts.maxStep);
        break;
    end

    if opts.debugMode
        o.averageAccuracy = 0.99 * o.averageAccuracy + 0.01 * ham.accuracy;
    end
    
    o.i = o.i + 1;
end
o.sampleTime = toc(t) - o.prepareTime;

opts.outputFunc('sample:end', 'Finished sampling. Now doing postprocessing.\n', o.i);
if (o.done == false)
    o.ess = effective_sample_size(o.samples);
    o.nSamples = min(o.ess);
    o.mixingTime = size(o.samples,2) / o.nSamples * o.iterPerRecord;
end

if ~opts.debugMode
    o.samples = P.T * o.samples + P.y;
end

%if isempty(P.df)
%	[o.pVal] = uniformtest(o.samples, problem, o.dim,  struct('toPlot', 1));
%end

% if isempty(P.df)
%     if opts.debugMode == true
%         s = o.ham.T * o.samples + o.ham.y;
%         [o.pVal] = uniformtest(s, problem, o.dim,  struct('toPlot', 1));
%     else
%         [o.pVal] = uniformtest(o.samples, problem, o.dim,  struct('toPlot', 1));
%     end
% end

opts.outputFunc('sample:end', 'See log.text for full output; samples in o.samples\n');
opts.outputFunc('sample:end', 'Mixing Time: %fs,  Effective Sample Size: %i\n', o.mixingTime, round(min(o.ess)));