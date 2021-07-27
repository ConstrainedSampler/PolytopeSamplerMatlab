function o = sample(problem, N, opts)
%Input: a structure P with the following fields
%  .Aineq
%  .bineq
%  .Aeq
%  .beq
%  .lb
%  .ub
%  .f
%  .df
%  .ddf
% describing a log-concave distribution given by
%   exp(-sum f_i(x_i))
%       over
%   {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}
% where f is given by a vector function of its 1-st, 2-nd, 3-rd
% derivative.
%
% Case 1: df is not defined
%   f(x) = 0.
% In this case, f, ddf must be empty.
%
% Case 2: df is a vector
%   f_i(x_i) = df_i x_i.
% In this case, f, ddf must be empty.
%
% Case 3: df is a function handle
%   f need to be defined as a function handle.
%   df need to be the derivative of f
%   ddf is optional. Providing ddf improves the mixing time.
% N - number of independent samples
% opts - sampling options
%
%Output:
% o - a structure containing the following properties:
%   samples - a cell of dim x N vectors containing each chain of samples
%   independent_samples - independent samples extracted (according to effective sample size)
%   prepareTime - time to pre-process the input (including find interior
%                 point, remove redundant constraints, reduce dimension etc.)
%   sampleTime - total sampling time in seconds (sum over all workers)
t = tic;

%% Initialize parameters and compiling if needed
if (nargin <= 2)
    opts = default_options();
else
    opts = setfield(default_options(), opts);
end
opts.startTime = t;

compile_solver(0); compile_solver(opts.simdLen);

%% Presolve
rng(opts.seed, 'simdTwister');
opts.seed = rng().Seed;

if ischar(opts.logging) || isstring(opts.logging) % logging for Polytope
    fid = fopen(opts.logging, 'a');
    opts.presolve.logFunc = @(tag, msg) fprintf(fid, '%s', msg);
elseif ~isempty(opts.logging)
    opts.presolve.logFunc = opts.logging;
else
    opts.presolve.logFunc = @(tag, msg) 0;
end

polytope = Polytope(problem, opts);

if ischar(opts.logging) || isstring(opts.logging)
    fclose(fid);
end

prepareTime = toc(t);


%% Check the trivial case
if polytope.n == 0
    warning('The domain consists only a single point.');
    o = struct;
    o.prepareTime = prepareTime;
    o.sampleTime = 0;
    o.problem = polytope;
    o.samples = polytope.center;
    return
end

%% Set up workers if nWorkers ~= 1
if opts.nWorkers ~= 1 && canUseParallelPool
    % create pool with size nWorkers
    p = gcp('nocreate');
    if isempty(p)
        if opts.nWorkers ~= 0
            p = parpool(opts.nWorkers);
        else
            p = parpool();
        end
    elseif opts.nWorkers ~= 0 && p.NumWorkers ~= opts.nWorkers
        delete(p);
        p = parpool(opts.nWorkers);
    end
    opts.nWorkers = p.NumWorkers;
    opts.N = N + opts.nRemoveInitialSamples * opts.nWorkers * opts.simdLen;
    
    spmd(opts.nWorkers)
        if opts.profiling
            mpiprofile on
		end
		
        rng(opts.seed + labindex, 'simdTwister');
        s = Sampler(polytope, opts);
        while s.terminate == 0
            s.step();
        end
        s.finalize();
        workerOutput = s.output;
        
        if opts.profiling
            mpiprofile viewer
        end
    end
    
    o = struct;
    o.workerOutput = cell(opts.nWorkers, 1);
    o.sampleTime = 0;
    for i = 1:opts.nWorkers
        o.workerOutput{i} = workerOutput{i};
        o.sampleTime = o.sampleTime + o.workerOutput{i}.sampleTime;
    end
    
    if ~opts.rawOutput
        o.chains = o.workerOutput{1}.chains;
        for i = 2:opts.nWorkers
            o.chains = [o.chains o.workerOutput{i}.chains];
            o.workerOutput{i}.chains = [];
        end
    end
else
    opts.N = N + opts.nRemoveInitialSamples * opts.simdLen;
	
    if opts.profiling
        profile on
    end
    
    s = Sampler(polytope, opts);
    while s.terminate == 0
        s.step();
    end
    s.finalize();
    o = s.output;
    
    if opts.profiling
        profile report
    end
end
o.problem = polytope;
o.opts = opts;

if ~opts.rawOutput
    o.ess = effective_sample_size(o.chains);
    
    y = [];
    for i = 1:numel(o.ess)
        chain_i = o.chains{i};
        ess_i = o.ess{i};
        N_i = size(chain_i,2);
        gap = ceil(N_i/ min(o.ess{i}, [], 'all'));
        for j = 1:size(ess_i,2)
            y_ij = chain_i(:, ceil(opts.nRemoveInitialSamples*gap:gap:N_i));
            y = [y y_ij];
        end
    end
    o.samples = y;
    o.summary = summary(o);
end

o.prepareTime = prepareTime;
