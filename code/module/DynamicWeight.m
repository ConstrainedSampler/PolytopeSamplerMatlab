classdef DynamicWeight < handle
    % Module for updating the weight in the log barrier.
    properties
       consecutiveBadStep = 0;
    end
    
    methods
        function o = DynamicWeight(s)
            extrHess = s.ham.barrier.extraHessian;
            s.ham.barrier = WeightedTwoSidedBarrier(s.ham.barrier.lb, s.ham.barrier.ub, s.problem.w, 2);
            s.ham.barrier.extraHessian = extrHess;
        end
        
        function o = step(o, s)
            bad_step = s.prob < 0.5 | s.ODEStep == s.opts.maxODEStep;
            o.consecutiveBadStep = bad_step .* o.consecutiveBadStep + bad_step;
            if (~s.freezed && ~all(s.accept))
                lsc = max(s.ham.lsc, [], 1);
                w = reshape(s.ham.barrier.w, size(lsc));
                if max(o.consecutiveBadStep) > 2
                   threshold = 4;
                else
                   threshold = 16;
                end
                
                idx = find(lsc > threshold * w);
                if ~isempty(idx)
                    s.ham.barrier.w(idx) = min(s.ham.barrier.w(idx) * threshold, 1);
                    s.ham.move(s.x, true); % Update the cache used in barrier
                    s.v = s.ham.resample(s.x, zeros(size(s.x)));
                    s.log('DynamicWeight:update_weight', 'The weight of %i coordinates are changed.\n', length(idx));
                end
            end
        end
    end
end