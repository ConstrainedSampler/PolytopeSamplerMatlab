classdef MixingTimeEstimator < handle
    % Module for estimate mixing time
    
    properties
        opts
        
        removedInitial = false
        sampleRate = 0
        sampleRateOutside = 0
        estNumSamples = 0
        estNumSamplesOutside = 0
        nextEstimateStep
    end
    
    methods
        function o = MixingTimeEstimator(s)
            o.opts = s.opts.MixingTimeEstimator;
            o.nextEstimateStep = o.opts.initialStep;
        end
        
        function o = step(o, s)
            if mean(s.nEffectiveStep) > o.nextEstimateStep
                ess = effective_sample_size(s.chains);
                ess = min(ess, [], 'all');
                
                if (o.removedInitial == false && ess > 2 * s.opts.nRemoveInitialSamples)
                    k = ceil(s.opts.nRemoveInitialSamples * (size(s.chains, 3) / ess));
                    s.chains = s.chains(:,:,k:end);
                    s.i = ceil(s.i * (1-k / size(s.chains, 3)));
                    o.removedInitial = true;
                    ess = effective_sample_size(s.chains);
                    ess = min(ess, [], 'all');
                end

                s.mixingTime = s.iterPerRecord * size(s.chains, 3) / ess;
                o.sampleRate = size(s.chains,1) / s.mixingTime;
                o.estNumSamples = s.i * o.sampleRate;
                s.share('sampleRate', o.sampleRate);
                s.share('estNumSamples', o.estNumSamples);
                o.update(s);
            end
        end
        
        function o = sync(o, s)
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
            o.update(s);
        end
        
        function o = update(o, s)
            s.sampleRate = o.sampleRate + o.sampleRateOutside;
            s.totalNumSamples = o.estNumSamples + o.estNumSamplesOutside;
            if o.estNumSamples > s.opts.freezeMCMCAfterSamples
                s.freezed = true;
            end
            
            if s.totalNumSamples > s.N && o.removedInitial
                s.share('estNumSamples', o.estNumSamples);
                s.terminate = 1;
                s.log('sample:end', '%i samples found.\n', s.totalNumSamples);
            else
                estimateEndingStep = s.N / s.sampleRate * (mean(s.nEffectiveStep) / s.i);
                o.nextEstimateStep = min(o.nextEstimateStep * o.opts.stepMultiplier, estimateEndingStep);
            end
        end
    end
end