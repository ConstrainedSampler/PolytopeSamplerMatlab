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
        freezed = false
        H1, H2
        ODEStep
        
        % Output related
        output = struct
        
        % Samples storaged for estimating number of samples
        samples = []
        mixingTime = NaN    % mixing time for this worker
        sampleRate = NaN    % sample rate achieved by all workers
        totalNumSamples = 0 % total number of samples by all workers
        iterPerRecord = 1
    end
    
    methods
        function log(o, tag, format, varargin)
            msg = sprintf('iter %i:', o.i);
            if nargin >= 4
                msg = [msg, sprintf(format, varargin{:})];
            else
                msg = [msg, sprintf(format)];
            end
            o.opts.logFunc(tag, msg, o);
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