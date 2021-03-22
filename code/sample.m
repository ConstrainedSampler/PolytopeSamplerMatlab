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
%   prepareTime - time to pre-process the input (including find interior
%                 point, remove redundant constraints, reduce dimension etc.)
%   sampleTime - total sampling time in seconds (sum over all workers)
t = tic;

%% Initialize parameters and compiling if needed
if (nargin <= 2)
    opts = default_options();
end
compile_solver(0); compile_solver(opts.simdLen); 
opts.startTime = t;
outputFormat = opts.outputFormat;

%% Presolve
rng(opts.seed, 'simdTwister');
polytope = Polytope(problem, opts);
prepareTime = toc(tic);

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
    nWorkers = p.NumWorkers;
    
    % generate random seeds
    seeds = randi(intmax('uint32'), nWorkers, 1, 'uint32');
    
    spmd(nWorkers)
        opts.seed = seeds(labindex);
        opts.labindex = labindex;
        workerOutput = sample_inner(problem, N, opts, polytope, nWorkers);
    end
    
    o = struct;
    o.workerOutput = cell(nWorkers, 1);
    o.sampleTime = 0;
    for i = 1:nWorkers
        o.workerOutput{i} = workerOutput{i};
        o.sampleTime = o.sampleTime + o.workerOutput{i}.sampleTime;
    end
    
    if ~strcmp(outputFormat,'raw')
        o.samples = o.workerOutput{1}.samples;
        for i = 2:nWorkers
            o.samples = [o.samples o.workerOutput{i}.samples];
            o.workerOutput{i}.samples = [];
        end
    end
else
    opts.seed = randi(intmax('uint32'), 'uint32');
    opts.labindex = 1;
    o = sample_inner(problem, N, opts, polytope, 1);
end

if ~strcmp(outputFormat,'raw')
    o.summary = summary(o.samples);
    ess = min(o.summary.ess);
    if iscell(o.samples)
        nRecord = 0;
        for i = 1:numel(o.samples)
            nRecord = nRecord + size(o.samples{i}, 2);
        end
    else
        nRecord = size(o.samples, 2);
    end
    o.mixingRecord = nRecord / ess;
end
o.prepareTime = prepareTime;
