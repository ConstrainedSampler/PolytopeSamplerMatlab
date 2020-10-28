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
%   mixingRecord - number of records per "independent" sample based on min(ess)
%   prepareTime - time to pre-process the input (including find analytic
%   center, remove redundant constraints, reduce dimension etc.)
%   sampleTime - total sampling time in seconds
t = tic;

%% Initialize parameters
if (nargin <= 2)
    opts = default_options();
end

s = Sampler;
s.problem = problem;
s.opts = opts;
r = rng(opts.seed, 'simdTwister');
s.seed = r.Seed;
s.N = N;
opts.module = unique([{'MixingTimeEstimator', 'SampleStorage'}, opts.module], 'stable');
for i = 1:length(opts.module)
    s.module{end+1} = feval(opts.module{i}, s);
end

%% Presolve
s.polytope = Polytope(problem, opts);
s.output.polytope = s.polytope;

%% Sample
s.ham = Hamiltonian(s.polytope, opts);
s.stepSize = opts.initalStepSize;
s.momentum = 1 - min(1, s.stepSize/opts.effectiveStepSize);
s.x = s.polytope.center;
s.v = s.ham.resample(s.x, zeros(size(s.x)));

for i = 1:length(opts.module)
    s.module{i}.initialize();
end

s.ham.opts.checkPrecision = true;
while true
    % v step
    s.v_ = s.ham.resample(s.x, s.v, s.momentum);
    
    % x step
    s.H1 = s.ham.H(s.x, s.v_);
    [s.x_, s.v_, s.ODEStep] = opts.odeMethod(s.x, s.v_, s.stepSize, s.ham, opts);
    s.feasible = s.ham.feasible(s.x_, s.v_);
    
    if ~s.feasible
        s.prob = 0;
    else
        s.H2 = s.ham.H(s.x_, -s.v_);
        s.prob = min(1, exp(s.H1-s.H2));
    end
    
    for i = 1:length(opts.module)
        s.module{i}.propose();
    end
    
    % rejection
    if rand >= s.prob
        s.log('sample:reject', 'rejected\n');
        s.v = -s.v; % flip the direction
        s.accept = false;
    else
        s.x = s.x_; s.v = s.v_;
        s.accept = true;
    end
    for i = 1:length(opts.module)
        s.module{i}.step();
    end
    
    if (toc(t) > opts.maxTime)
        s.log('sample:end', 'opts.maxTime (%f sec) reached.\n', opts.maxTime);
        s.terminate = 2;
    end
    
    if s.i >= opts.maxStep
        s.log('sample:end', 'opts.maxStep (%i step) reached.\n', opts.maxStep);
        s.terminate = 2;
    end
    
    if s.terminate > 0
        break;
    end
    
    s.i = s.i + 1;
    % check solver precision if not accepted or too many steps
    s.ham.opts.checkPrecision = (s.accept == false) || (s.ODEStep >= opts.maxODEStep);
end

for i = 1:length(opts.module)
    s.module{i}.finalize();
end

o = s.output;
% if (o.done == false)
%     o.ess = effective_sample_size(o.samples);
%     o.nSamples = min(o.ess);
%     o.mixingIter = size(o.samples,2) / o.nSamples * o.iterPerRecord;
% end
% o.mixingRecord = size(o.samples,2) / o.nSamples;
% 
% opts.outputFunc('sample:end', 'See log.text for full output; samples in o.samples\n');
% opts.outputFunc('sample:end', 'Mixing Time (in iter): %f,  Effective Sample Size: %i\n', o.mixingIter, round(min(o.ess)));