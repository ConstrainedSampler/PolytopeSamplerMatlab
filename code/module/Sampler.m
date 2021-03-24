classdef Sampler < dynamicprops
    properties
        % Problem formulation
        problem, N, opts % inputs
        polytope
        seed % random seed used in the sampler
        
        % System related
        module = {}
        i = 1
        terminate = 0       % 0 = not terminated, 
                            % 1 = terminated with enough samples,
                            % 2 = timeout / too many steps,
                            % 3 = error
        nWorkers
        labindex
        shareLog = struct()   % variables shared to different workers
        logUpdated = false    % if share variable is updated
        lastBroadcast         % time last broadcast
        shared = {}           % variables shared from different workers
        
        
        % HMC parameters
        ham
        stepSize
        momentum
        
        % HMC state
        x, v   % current point
        x_, v_ % proposed point
        prob, accept, feasible
        acceptedStep = 0;
        freezed = false
        H1, H2
        ODEStep
        
        % Output related
        output = struct
        fid
        
        % Samples storaged for estimating number of samples
        samples = []
        mixingTime = NaN    % mixing time for this worker
        sampleRate = NaN    % sample rate achieved by all workers
        totalNumSamples = 0 % total number of samples by all workers
        iterPerRecord = 1
    end
    
    methods
        function o = Sampler(problem, polytope, nWorkers, opts)
            o.problem = problem;
            o.polytope = polytope;
            o.opts = opts;
            o.nWorkers = nWorkers;
            o.seed = opts.seed;
            o.N = opts.N;
            rng(opts.seed, 'simdTwister');
            
            if (o.nWorkers > 1)
                o.labindex = labindex;
                o.lastBroadcast = tic;
                for i = 1:nWorkers
                    o.shared{i} = struct;
                end
            else
                o.labindex = 1;
            end
            
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
            
            opts.module = unique([{'MixingTimeEstimator', 'SampleStorage'}, opts.module], 'stable');
            for i = 1:length(opts.module)
                o.module{end+1} = feval(opts.module{i}, o);
            end
        end
        
        function finalize(o)
            for j = 1:length(o.opts.module)
                o.module{j}.finalize();
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
            
            if isa(o.opts.logging,'function_handle')
                o.opts.logging(tag, msg, o);
            else
                fprintf(o.fid, '%s', msg);
            end
        end
        
        function share(o, name, var)
            if (o.nWorkers > 1)
                o.logUpdated = true;
                o.shareLog.(name) = var;
                o.shared{o.labindex}.(name) = var;
            end
        end
    end
end