% TODO: Understand why I need to multiply by n.
% Maybe we check TV ball to see what is the smallest I can do

classdef DynamicRegularizer
    properties
        sampler
        prerequisites = {}
        
        bound
    end
    
    methods
        function o = DynamicRegularizer(sampler)
            o.sampler = sampler;
            o.setBound(max(abs(sampler.x), 1));
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