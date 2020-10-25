classdef DebugLogger < handle
    properties
        sampler
        startTime
    end
    
    methods
        function o = DebugLogger(sampler)
            o.sampler = sampler;
            o.startTime = tic();
            sampler.output.acceptedStep = 0;
            sampler.output.averageLinearSystemAccuracy = 0;
        end
        
        function o = initialize(o)
            o.sampler.output.prepareTime = toc(o.startTime);
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            s.output.acceptedStep = s.output.acceptedStep + s.accept;
            s.output.averageLinearSystemAccuracy = 0.9 * s.output.averageLinearSystemAccuracy + 0.1 * o.sampler.ham.accuracy;
        end
        
        function o = finalize(o)
            o.sampler.output.sampleTime = toc(o.startTime) - o.sampler.output.prepareTime;
            o.sampler.output.sampler = o.sampler;
            
            if o.sampler.terminate == 3
                save('dump.mat');
            end
        end
    end
end