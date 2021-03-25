function o = sample(problem, N, opts)
%Input:
% problem - a struct containing the constraints defining the polytope
%
% N - number of indepdent samples
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
end
opts.startTime = t;
opts.N = N;

compile_solver(0); compile_solver(opts.simdLen);

%% Presolve
rng(opts.seed, 'simdTwister');

if ischar(opts.logging) || isstring(opts.logging) % logging for Polytope
    fid = fopen(opts.logging, 'a');
    opts.logFunc = @(tag, msg) fprintf(fid, '%s', msg);
elseif ~isempty(opts.logging)
    opts.logFunc = opts.logging;
else
    opts.logFunc = @(tag, msg) 0;
end

polytope = Polytope(problem, opts);

if ischar(opts.logging) || isstring(opts.logging)
    fclose(fid);
end

prepareTime = toc(t);

%% Set up workers if nWorkers ~= 1
if opts.nWorkers ~= 1
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
    
    % generate random seeds
    seeds = randi(intmax('uint32'), nWorkers, 1, 'uint32');
    
    spmd(opts.nWorkers)
        if opts.profiling
            mpiprofile on
        end
        
        opts_local = opts;
        opts_local.seed = seeds(labindex);
        workerOutput = sample_inner(polytope, opts_local);
        
        if opts.profiling
            mpiprofile viewer
        end
    end
    
    o = struct;
    o.workerOutput = cell(nWorkers, 1);
    o.sampleTime = 0;
    for i = 1:nWorkers
        o.workerOutput{i} = workerOutput{i};
        o.sampleTime = o.sampleTime + o.workerOutput{i}.sampleTime;
    end
    
    if ~opts.rawOutput
        o.samples = o.workerOutput{1}.chains;
        for i = 2:nWorkers
            o.samples = [o.samples o.workerOutput{i}.chains];
            o.workerOutput{i}.chains = [];
        end
    end
else
    opts.seed = randi(intmax('uint32'), 'uint32');
    opts.nWorkers = 1;
    if opts.profiling
        profile on
    end
    
    o = sample_inner(polytope, opts);
    
    if opts.profiling
        profile report
    end
end

if ~opts.rawOutput
    o.ess = effective_sample_size(o.samples);
    o.summary = summary(o);
    
    y = [];
    for i = 1:numel(o.ess)
        samples_i = o.samples{i};
        ess_i = o.ess{i};
        N_i = size(samples_i,2);
        gap = ceil(N_i/ min(o.ess{i}, [], 'all'));
        for j = 1:size(ess_i,2)
            y_ij = samples_i(:, ceil(opts.nRemoveInitialSamples*gap:gap:N_i));
            y = [y y_ij];
        end
    end
    o.independent_samples = y;
end

o.prepareTime = prepareTime;
o.polytope = polytope;
