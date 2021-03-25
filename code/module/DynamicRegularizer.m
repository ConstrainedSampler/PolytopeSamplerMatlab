classdef DynamicRegularizer < handle
    % Module for updating the extra term we add to the barrier
    % This is nesscary for any polytope with free variables
    properties
        bound = 1
    end
    
    methods
        function o = DynamicRegularizer(s)
            s.ham.barrier.extraHessian = 1;
            o.set_bound(max(abs(s.x)), s);
        end
        
        function o = propose(o, s)
            o.set_bound(max(abs(s.x)), s);
        end
        
        function o = set_bound(o, bound, s)
            o.bound = max(bound, o.bound);
            if (~s.freezed)
                idx = find(1./(bound.*bound) < s.ham.n * s.ham.barrier.extraHessian);
                if ~isempty(idx)
                    s.ham.barrier.extraHessian = 0.2./(s.ham.n * o.bound.*o.bound);
                    s.v = s.ham.resample(s.x, zeros(size(s.x)));
                    s.log('DynamicRegularizer:set_bound', 'The bound of %i coordinates are changed.\n', length(idx));
                end
            end
        end
    end
end