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
        lastBroadcast           % time last broadcast
        everyOneElse            % array [nWorkers] \setminus  labindex
        
        % System properties
        module = {}
        startTime               % starting time for sampling
        terminate = 0           % 0 = not terminated,
                                % 1 = terminated with enough samples,
                                % 2 = timeout / too many steps,
                                % 3 = error
                                % <0 = terminated by other samples
        
        
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
        
        % Statiatics
        accumulatedMomentum = 0 % accumulated momentum
        nEffectiveStep = 0      % number of effective steps
        acceptedStep = 0
    end
    
    methods
        function o = Sampler(problem, opts)
            o.problem = problem;
            o.opts = opts;
            o.nWorkers = opts.nWorkers;
            o.N = opts.N;
            o.startTime = tic();
            
            % Initialize parallel properties
            if (o.nWorkers > 1)
                o.labindex = labindex;
                for j = 1:o.nWorkers
                    o.shared{j} = struct();
                end
                
                o.everyOneElse = 1:o.nWorkers;
                o.everyOneElse(o.labindex) = [];
                o.lastBroadcast = tic();
            end
            rng(bitxor(opts.seed, uint32(o.labindex)), 'simdTwister');
            
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
            o.x = ones(opts.simdLen, 1) * problem.center';
            
            % Initialize modules
            opts.module = unique([{'MixingTimeEstimator', 'MemoryStorage'}, opts.module], 'stable');
            for j = 1:length(opts.module)
                o.module{end+1} = feval(opts.module{j}, o);
            end
            
            % Initialize the point
            o.v = o.ham.resample(o.x, zeros(size(o.x)));
        end
        
        function step(o)
            % v step
            o.v_ = o.ham.resample(o.x, o.v, o.momentum);
            
            % x step
            o.H1 = o.ham.H(o.x, o.v_);
            [o.x_, o.v_, o.ODEStep] = o.opts.odeMethod(o.x, o.v_, o.stepSize, o.ham, o.opts);
            o.feasible = o.ham.feasible(o.x_, o.v_);
            o.H2 = o.ham.H(o.x_, -o.v_);
            o.prob = min(1, exp(o.H1-o.H2)) .* o.feasible;

            % rejection
            o.accept = (rand(o.opts.simdLen, 1) < o.prob);
            o.x = blendv(o.x, o.x_, o.accept);
            o.v = blendv(-o.v, o.v_, o.accept);

            if (~all(o.accept))
                o.log('sample:reject', 'rejected chain ' + join(string(num2str(find(~o.accept))),',') + '\n');
            end
            
            o.acceptedStep = o.acceptedStep + o.prob;
            o.accumulatedMomentum = mean(o.prob) * o.momentum * o.accumulatedMomentum + o.stepSize;
            o.nEffectiveStep = o.nEffectiveStep + o.stepSize * o.accumulatedMomentum * o.accept;
            for j = 1:length(o.opts.module)
                o.module{j}.step(o);
            end
            
            if (toc(o.opts.startTime) > o.opts.maxTime)
                o.log('sample:end', 'opts.maxTime (%f sec) reached.\n', o.opts.maxTime);
                o.terminate = 2;
                return;
            end

            if o.i >= o.opts.maxStep
                o.log('sample:end', 'opts.maxStep (%i step) reached.\n', o.opts.maxStep);
                o.terminate = 2;
                return;
            end

            if (o.nWorkers > 1)
                o.recieve();
                o.broadcast();
            end
            
            o.i = o.i + 1;
        end
        
        function finalize(o)
            for j = 1:length(o.opts.module)
                if ismethod(o.module{j}, 'finalize')
                    o.module{j}.finalize(o);
                end
            end
            
            if isstring(o.opts.logging) || ischar(o.opts.logging)
                fclose(o.fid);
            end
            
            if (o.nWorkers > 1 && o.terminate > 0)
                o.shareCache.terminate = o.terminate;
                labSend(o.shareCache, o.everyOneElse);
            end
            
            if (o.nWorkers > 1) % avoid the warning (An incoming message was discarded)
                labBarrier;
                while (labProbe)
                    labReceive;
                end
            end
            
            o.output.sampleTime = toc(o.startTime);
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
        
        function broadcast(o)
            if (toc(o.lastBroadcast) > o.opts.broadcastInterval && ~isempty(fieldnames(o.shareCache)))
                labSend(o.shareCache, o.everyOneElse);
                o.shareCache = struct();
                o.lastBroadcast = tic();
            end
        end
        
        function recieve(o)
            updated = false;
            while (labProbe)
                [data, idx, ~] = labReceive;
                
                f = fieldnames(data);
                for j = 1:length(f)
                    o.shared{idx}.(f{j}) = data.(f{j});
                end
                
                if isfield(data, 'terminate')
                    o.terminate = -data.terminate;
                end
                updated = true;
            end
            
            if (updated)
                for j = 1:length(o.module)
                    if ismethod(o.module{j}, 'sync')
                        o.module{j}.sync(o);
                    end
                end
            end
        end
    end
end