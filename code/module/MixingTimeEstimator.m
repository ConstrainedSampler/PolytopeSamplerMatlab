classdef MixingTimeEstimator < handle
    properties
        sampler
        opts
        
        nextEstimateIter
        ess
    end
    
    methods
        function o = MixingTimeEstimator(sampler)
            o.sampler = sampler;
            o.opts = sampler.opts.MixingTimeEstimator;
            o.nextEstimateIter = o.opts.startIter;
        end
        
        function o = initialize(o)
            
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            if s.i > o.nextEstimateIter
                o.ess = effective_sample_size(s.samples);
                s.mixingTime = s.iterPerRecord * size(s.samples,2) / min(o.ess);
                lastEstNumSamples = s.i / s.mixingTime;
                
                if lastEstNumSamples > s.N
                    s.terminate = 1;
                    s.log('sample:end', '%i samples found.\n', lastEstNumSamples);
                end
                
                estimateEndingIter = s.i + (s.N - lastEstNumSamples) * s.mixingTime;
                o.nextEstimateIter = min(o.nextEstimateIter * o.opts.iterMultiplier, estimateEndingIter);
            end
        end
        
        function o = finalize(o)
            
        end
    end
end