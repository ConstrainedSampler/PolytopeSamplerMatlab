function o = sample(problem, N, opts)
%Input: a structure problem with the following fields
%  .Aineq
%  .bineq
%  .Aeq
%  .beq
%  .lb
%  .ub
%  .f
%  .df
%  .ddf
% describing a logconcave distribution given by
%   exp(-sum f_i(x_i))
%       over
%   {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}
% where f is given by a vector function of its first (df) and second (ddf)
% derivatives.
%
% Case 1: df is not defined
%   f(x) = 0.
% In this case, f, ddf must be empty. This is the uniform distribution. In
% this case the feasible region must be bounded.
%
% Case 2: df is a vector
%   f_i(x_i) = df_i x_i.
% In this case, f, ddf must be empty.
%
% Case 3: df is a function handle
%   f need to be defined as a function handle.
%   df need to be the derivative of f
%   ddf is optional. Providing ddf could improve the mixing time.
%
% N - number of independent samples 
%
% opts - sampling options
%
%Output:
% o - a structure containing the following properties:
%   samples - a cell of dim x N vectors containing each chain of samples
%   prepareTime - time to pre-process the input (including find interior
%                 point, remove redundant constraints, reduce dimension)
%   sampleTime - total sampling time in seconds (sum over all workers)
t = tic;

%% Initialize parameters and compiling if needed
if (nargin <= 2)
    opts = default_options();
else
    opts = Setfield(default_options(), opts);
end
opts.startTime = t;

compile_solver(0); compile_solver(opts.simdLen);

%% Presolve
if isempty(opts.seed)
    opts.seed = randi(2^31);
end

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
if opts.nWorkers ~= 1 && ~isempty(ver('parallel'))
    % create pool with size nWorkers
    p = gcp('nocreate');
    if isempty(p)
        if opts.nWorkers ~= Inf
            p = parpool(opts.nWorkers);
        else
            p = parpool();
        end
    elseif opts.nWorkers ~= Inf && p.NumWorkers ~= opts.nWorkers
        delete(p);
        p = parpool(opts.nWorkers);
    end
    opts.nWorkers = p.NumWorkers;
    opts.N = N;
    
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
    o.prepareTime = prepareTime;
    for i = 1:opts.nWorkers
        o.workerOutput{i} = workerOutput{i};
        o.sampleTime = o.sampleTime + o.workerOutput{i}.sampleTime;
        o.prepareTime = o.prepareTime + o.workerOutput{i}.prepareTime;
    end
    
    if ~opts.rawOutput
        o.chains = o.workerOutput{1}.chains;
        for i = 2:opts.nWorkers
            o.chains = [o.chains o.workerOutput{i}.chains];
            o.workerOutput{i}.chains = [];
        end
    end
else
    opts.N = N;
	
    if opts.profiling
        profile on
    end
    
    rng(opts.seed, 'simdTwister');
    s = Sampler(polytope, opts);
    while s.terminate == 0
        s.step();
    end
    s.finalize();
    o = s.output;
    o.prepareTime = o.prepareTime + prepareTime;
    
    if opts.profiling
        profile report
    end
end
o.problem = polytope;
o.opts = opts;

if ~opts.rawOutput
    o.samples = o.chains;
    o = rmfield(o, 'chains');
    o.summary = summary(o);
end
