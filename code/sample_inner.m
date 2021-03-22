function o = sample_inner(problem, N, opts, polytope, nWorkers)    %% Set up Sampler
    startTime = tic;
    s = Sampler;
    s.problem = problem;
    s.opts = opts;
    s.seed = opts.seed;
    s.nWorkers = nWorkers;
    if (s.nWorkers > 1)
        s.labindex = labindex;
        everyOneElse = 1:nWorkers;
        everyOneElse(labindex) = [];
        s.lastBroadcast = tic;
        for i = 1:opts.nWorkers
            s.shared{i} = struct;
        end
    else
        s.labindex = 1;
    end
    rng(opts.seed, 'simdTwister');
    s.N = N;
    s.polytope = polytope;
    s.output.polytope = s.polytope;
    opts.module = unique([{'MixingTimeEstimator', 'SampleStorage'}, opts.module], 'stable');
    for i = 1:length(opts.module)
        s.module{end+1} = feval(opts.module{i}, s);
    end
    
    %% Sample
    s.ham = Hamiltonian(s.polytope, opts);
    s.stepSize = opts.initalStepSize;
    s.momentum = 1 - min(1, s.stepSize/opts.effectiveStepSize);
    s.x = ones(opts.simdLen, 1) * s.polytope.center';
    s.v = s.ham.resample(s.x, zeros(size(s.x)));

    for i = 1:length(opts.module)
        s.module{i}.initialize();
    end

    while true
        % v step
        s.v_ = s.ham.resample(s.x, s.v, s.momentum);

        % x step
        s.H1 = s.ham.H(s.x, s.v_);
        [s.x_, s.v_, s.ODEStep] = opts.odeMethod(s.x, s.v_, s.stepSize, s.ham, opts);
        s.feasible = s.ham.feasible(s.x_, s.v_);
        for i = 1:length(opts.module)
            s.module{i}.propose();
        end

        s.H2 = s.ham.H(s.x_, -s.v_);
        s.prob = min(1, exp(s.H1-s.H2)) .* s.feasible;

        % rejection
        s.accept = (rand(opts.simdLen, 1) < s.prob);
        s.x = blendv(s.x, s.x_, s.accept);
        s.v = blendv(-s.v, s.v_, s.accept);

        if (~all(s.accept))
            s.log('sample:reject', 'rejected chain ' + join(string(num2str(find(~s.accept))),',') + '\n');
        end

        for i = 1:length(opts.module)
            s.module{i}.step();
        end

        if (toc(opts.startTime) > opts.maxTime)
            s.log('sample:end', 'opts.maxTime (%f sec) reached.\n', opts.maxTime);
            s.terminate = 2;
        end

        if s.i >= opts.maxStep
            s.log('sample:end', 'opts.maxStep (%i step) reached.\n', opts.maxStep);
            s.terminate = 2;
        end

        if s.terminate > 0
            if (s.nWorkers > 1)
                s.shareLog.terminate = s.terminate;
                labSend(s.shareLog, everyOneElse);
            end
            break;
        end
        
        if (s.nWorkers > 1)
            % receive
            updated = false; terminateReceived = false;
            while (labProbe)
                [data,idx,~] = labReceive;
                f = fieldnames(data);
                for i = 1:length(f)
                    s.shared{idx}.(f{i}) = data.(f{i});
                end
                
                if isfield(data, 'terminate')
                    s.terminate = data.terminate;
                    terminateReceived = true;
                end
                updated = true;
            end
            
            if (updated)
                for i = 1:length(opts.module)
                    s.module{i}.sync();
                end
            end
            
            if terminateReceived % no broadcast if terminated by others
                break;
            end

            % broadcast
            if (toc(s.lastBroadcast) > opts.broadcastInterval && s.logUpdated)
                labSend(s.shareLog, everyOneElse);
                s.logUpdated = false;
                s.shareLog = struct;
            end
        end
        s.i = s.i + 1;
    end

    for i = 1:length(opts.module)
        s.module{i}.finalize();
    end
    
    o = s.output;
    o.sampleTime = toc(startTime);
    
    
    if (s.nWorkers > 1) % just to avoid warning (An incoming message was discarded)
        labBarrier;
        while (labProbe)
            labReceive;
        end
    end
end
