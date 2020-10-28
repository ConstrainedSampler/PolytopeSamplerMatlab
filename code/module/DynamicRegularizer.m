% TODO: Understand why I need to multiply by n.
% Maybe we check TV ball to see what is the smallest I can do

classdef DynamicRegularizer < handle
    properties
        sampler
        bound
    end
    
    methods
        function o = DynamicRegularizer(sampler)
            o.sampler = sampler;
        end
        
        function o = initialize(o)
            o.sampler.ham.barrier.extraHessian = 1;
            o.setBound(max(abs(o.sampler.x), 1));
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
                    s.ham.barrier.extraHessian = 0.1./(s.ham.n * o.bound.*o.bound);
                    s.v = s.ham.resample(s.x, zeros(size(s.x)));
                    s.log('DynamicWeight:setBound', 'The bound of %i coordinates are changed significantly.\n', length(idx));
                end
            end
        end
    end
end