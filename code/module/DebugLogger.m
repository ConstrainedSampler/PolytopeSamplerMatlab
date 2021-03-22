classdef DebugLogger < handle
    properties
        sampler
    end
    
    methods
        function o = DebugLogger(sampler)
            o.sampler = sampler;
            sampler.output.acceptedStep = 0;
            sampler.output.totalStep = 0;
            sampler.output.averageAccuracy = 0;
        end
        
        function o = initialize(o)
            
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            s.output.acceptedStep = s.output.acceptedStep + s.accept;
            s.output.totalStep = s.output.totalStep + 1;
            s.output.averageAccuracy = s.output.averageAccuracy + o.sampler.ham.solver.accuracy;
            if (any(o.sampler.ham.solver.accuracy > o.sampler.opts.solverThreshold))
                s.log('Solver:inaccurate', 'low double accuracy %s.\n', num2str(o.sampler.ham.solver.accuracy'));
            end
        end
        
        function o = sync(o)
            
        end
        
        function o = finalize(o)
            s = o.sampler;
            s.output.numCholesky = s.ham.solver.getDecomposeCount();
            s.output.averageAccuracy = s.output.averageAccuracy / s.output.totalStep;
            
            s.output.sampler = o.sampler; %TODO
            if s.nWorkers == 1
                s.output.sampler = o.sampler;
                if s.terminate == 3
                    save('dump_ignore.mat');
                end
            end
        end
    end
end