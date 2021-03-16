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
            sampler.output.totalStep = 0;
        end
        
        function o = initialize(o)
            o.sampler.output.prepareTime = toc(o.startTime);
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            s.output.acceptedStep = s.output.acceptedStep + s.accept;
            s.output.totalStep = s.output.totalStep + 1;
        end
        
        function o = finalize(o)
            s = o.sampler;
            s.output.nChol = s.ham.solver.getDecomposeCount();
            s.output.sampleTime = toc(o.startTime) - o.sampler.output.prepareTime;
            s.output.sampler = o.sampler;
            
            if s.terminate == 3
                save('dump.mat');
            end
        end
    end
end