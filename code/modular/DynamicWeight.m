% TODO: Understand why this fails
% if opts.weightedBarrier
%     o.averageLSC = zeros(size(o.x));
% end
%     if opts.weightedBarrier && o.nSamples < opts.freezeMCMCAfterSamples
%         o.averageLSC = 0.9 * o.averageLSC + 0.1 * ham.last_lsc;
%         idx = find(o.averageLSC > (0.5 / o.stepSize) * ham.barrier.weight);
%         if ~isempty(idx)
%             opts.outputFunc('sample:weight_change', 'iter %i: changed the weight of %i coordinates.\n', o.i, length(idx));
%             ham.barrier.weight(idx) = min(1, ham.barrier.weight(idx) * 2);
%             o.averageLSC = zeros(size(o.x));
%         end
%     end

% 
% classdef DynamicWeight
%     properties
%         sampler
%         opts
%         prerequisites = {}
%     end
%     
%     methods
%         function o = DynamicWeight(sampler)
%             o.sampler = sampler;
%             o.opts = sampler.opts.DynamicWeight;
%         end
%         
%         function o = propose(o)
%             
%         end
%         
%         function o = step(o)
%             
%         end
%         
%         function o = finalize(o)
%             
%         end
%     end
% end