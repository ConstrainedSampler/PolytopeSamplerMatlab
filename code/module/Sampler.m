classdef Sampler < dynamicprops
    properties
        % Input
        problem                 % Polytope object
        N                       % Number of samples needed
        opts                    % Options
        
        % Parallel properties
        nWorkers = 1            % Number of workers
        labindex = 1            % The index of this sampler
        shared = {}             % shared{i} is the variables shared by worker i
        shareCache = struct()   % variables that will be shared to different workers
        lastBroadcast = tic()   % time last broadcast
        
        % System properties
        module = {}
        startTime = tic()       % starting time for sampling
        terminate = 0           % 0 = not terminated,
        % 1 = terminated with enough samples,
        % 2 = timeout / too many steps,
        % 3 = error
        
        
        % HMC properties
        ham                     % Hamiltonian object
        stepSize                % step size in (0, inf)
        momentum                % momentum coefficient in [0, 1)
        freezed = false         % if markov chain is freezed (i.e. no parameter can changed)
        
        % HMC state
        i = 1                   % the iteration number
        x, v                    % current point
        x_, v_                  % proposed point
        prob                    % acceptance probability for the proposed point
        accept                  % if we accept the proposed point
        feasible                % if the proposed step is feasible
        ODEStep                 % the number of ODE steps we take for this proposed point
        H1, H2                  % Hamiltonian energy before and after
        
        % Output properties
        output = struct         % the return object
        fid                     % fileID for logging
        
        % MixingTime properties
        chains = []             % chains used to estimate mixing time,
        % depending on storage, it may not contains every samples
        mixingTime = NaN        % mixing time for this worker
        sampleRate = NaN        % sample rate achieved by all workers
        totalNumSamples = 0     % total number of samples by all workers
        iterPerRecord = 1       % iteration per record in chains
        acceptedStep = 0        % number of accepted steps
    end
    
    methods
        function o = Sampler(problem, opts)
            o.problem = problem;
            o.opts = opts;
            o.nWorkers = opts.nWorkers;
            o.N = opts.N;
            rng(opts.seed, 'simdTwister');
            
            % Initialize parallel properties
            if (o.nWorkers > 1)
                o.labindex = labindex;
                for i = 1:o.nWorkers
                    o.shared{i} = struct();
                end
            end
            
            % Initialize logging properties
            if isstring(o.opts.logging) || ischar(o.opts.logging)
                file_name = string(o.opts.logging);
                if (o.nWorkers > 1)
                    if endsWith(file_name, '.log')
                        file_name = extractBetween(file_name, 1, strlength(file_name)-4);
                        file_name = sprintf("%s_%i.log", file_name, labindex);
                    else
                        file_name = sprintf("%s%i", file_name, labindex);
                    end
                end
                o.fid = fopen(file_name, 'a');
            end
            
            % Initialize HMC properties
            o.ham = Hamiltonian(problem, opts);
            o.stepSize = opts.initalStepSize;
            o.momentum = 1 - min(1, o.stepSize/opts.effectiveStepSize);
            
            % Initialize modules
            opts.module = unique([{'MixingTimeEstimator', 'MemoryStorage'}, opts.module], 'stable');
            for i = 1:length(opts.module)
                o.module{end+1} = feval(opts.module{i}, o);
            end
            
            % Initialize the point
            o.x = ones(opts.simdLen, 1) * problem.center';
            o.v = o.ham.resample(o.x, zeros(size(o.x)));
        end
        
        function finalize(o)
            for j = 1:length(o.opts.module)
                if ismethod(o.module{j}, 'finalize')
                    o.module{j}.finalize();
                end
            end
            
            if isstring(o.opts.logging) || ischar(o.opts.logging)
                fclose(o.fid);
            end
        end
        
        function log(o, tag, format, varargin)
            if isempty(o.opts.logging)
                return;
            end
            
            msg = sprintf('iter %i:', o.i);
            if nargin >= 4
                msg = [msg, sprintf(format, varargin{:})];
            else
                msg = [msg, sprintf(format)];
            end
            
            if isa(o.opts.logging, 'function_handle')
                o.opts.logging(tag, msg, o);
            else
                fprintf(o.fid, '%s', msg);
            end
        end
        
        function share(o, name, var)
            if (o.nWorkers > 1)
                o.shareCache.(name) = var;
                o.shared{o.labindex}.(name) = var;
            end
        end
    end
end