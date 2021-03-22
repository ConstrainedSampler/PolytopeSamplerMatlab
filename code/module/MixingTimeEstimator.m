classdef MixingTimeEstimator < handle
    properties
        sampler
        opts
        
        sampleRate = 0
        sampleRateOutside = 0
        estNumSamples = 0
        estNumSamplesOutside = 0
        nextEstimateIter
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
                ess = effective_sample_size(s.samples);
                s.mixingTime = s.iterPerRecord * size(s.samples,3) / min(ess, [], 'all');
                o.sampleRate = size(s.samples,1) / s.mixingTime;
                o.estNumSamples = (s.i / s.mixingTime) * size(s.samples,1);
                s.share('sampleRate', o.sampleRate);
                s.share('estNumSamples', o.estNumSamples);
                o.check();
            end
        end
        
        function o = sync(o)
            s = o.sampler;
            o.estNumSamplesOutside = 0;
            o.sampleRateOutside = 0;
            for i = 1:s.nWorkers
                if i ~= s.labindex
                    if isfield(s.shared{i}, 'estNumSamples')
                        o.estNumSamplesOutside = o.estNumSamplesOutside + s.shared{i}.estNumSamples;
                    end
                    
                    if isfield(s.shared{i}, 'sampleRate')
                        o.sampleRateOutside = o.sampleRateOutside + s.shared{i}.sampleRate;
                    end
                end
            end
            o.check();
        end
        
        function o = check(o)
            s = o.sampler;
            s.totalNumSamples = o.estNumSamples + o.estNumSamplesOutside;
            if s.totalNumSamples > s.N
                s.share('estNumSamples', o.estNumSamples);
                %disp(labindex)
                s.terminate = 1;
                s.log('sample:end', '%i samples found.\n', s.totalNumSamples);
            end
            s.sampleRate = o.sampleRate + o.sampleRateOutside;
            estimateEndingIter = s.N / s.sampleRate;
            o.nextEstimateIter = min(o.nextEstimateIter * o.opts.iterMultiplier, estimateEndingIter);
        end
        
        function o = finalize(o)
            
        end
    end
end