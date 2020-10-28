classdef Sampler < dynamicprops
    properties
        % Problem formulation
        problem, N, opts % inputs
        polytope
        seed % random seed used in the sampler
        
        % System related
        module = {}
        i = 1
        terminate = 0	% 0 = not terminated, 
                        % 1 = terminated with enough samples,
                        % 2 = timeout / too many steps,
                        % 3 = error
        
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
        mixingTime = NaN
        iterPerRecord = 1
        
        % Logging related (Move to the class)
        startTime
        prepareTime
        sampleTime
        acceptedStep
    end
    
    methods
        function log(o, tag, format, varargin)
            msg = sprintf('iter %i:', o.i);
            if nargin == 4
                msg = [msg, sprintf(format, varargin{:})];
            else
                msg = [msg, sprintf(format)];
            end
            o.opts.logFunc(tag, msg);
        end
    end
end