%     if (o.i >= o.checkESSIter && o.acceptedSteps > opts.minAcceptedSteps)
%         len = size(o.samples,2);
%         o.ess = effective_sample_size(o.samples);
%         o.nSamples = min(o.ess);
%         o.mixingTime = o.i / o.nSamples;
%         
%         % thinning
%         if opts.recordsPerIndependentSample < Inf
%             h = max((len / o.nSamples) / opts.recordsPerIndependentSample, 1);
%             idx = ceil(h*(1:floor(len / h))) + (len - ceil(h * floor(len / h)));
%             o.samples = o.samples(:, idx);
%             o.iterPerRecord = ceil(o.mixingTime / (opts.recordsPerIndependentSample * 2));
%             if o.nSamples < N/2
%                 o.iterPerRecord = ceil(o.iterPerRecord / 2);
%             end
%         end
%         
%         if o.nSamples > N
%             o.done = true;
%             break;
%         end
%         
%         % check ess after "freezeMCMCAfterSamples" samples and N/2 samples
%         N_target = N;
%         if o.nSamples < N/2
%             N_target = min(N_target, N/2);
%         end
%         
%         if o.nSamples < opts.freezeMCMCAfterSamples
%             N_target = min(N_target, opts.freezeMCMCAfterSamples);
%         end
%         
%         o.checkESSIter = o.i + (N_target - o.nSamples) * o.mixingTime;
%     end

% I need to figure out how to update the ess....

classdef MixingTimeEstimator
    properties
        sampler
        prerequisites = {}
        
        bound
    end
    
    methods
        function o = DynamicRegularizer(sampler)
            o.sampler = sampler;
        end
        
        function o = propose(o)
            o.setBound(max(o.bound, abs(o.sampler.x)));
        end
        
        function o = step(o)
            
        end
        
        function o = finalize(o)
            
        end
        
        function o = setBound(o, bound)
            o.bound = bound;
            if (~o.sampler.freezed)
                s = o.sampler;
                idx = find(1./(bound.*bound) < s.ham.n * s.ham.barrier.extraHessian);
                if ~isempty(idx)
                    s.ham.barrier.extraHessian = 0.25./(s.ham.n * o.bound.*o.bound);
                    s.v = s.ham.resample(s.x, zeros(size(s.x)));
                    s.outputFunc('DynamicWeight:setBound', 'iter %i: the bound of %i coordinates are changed significantly.\n', s.i, length(idx));
                end
            end
        end
    end
end