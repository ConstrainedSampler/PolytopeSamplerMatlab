classdef DynamicRegularizer < handle
    % Module for updating the extra term we add to the barrier
    % This is nesscary for any polytope with free variables
    properties
        bound = 1
    end
    
    methods
        function o = DynamicRegularizer(s)
            s.ham.barrier.extraHessian = 1;
        end
        
        function o = step(o, s)
            o.bound = max(max(abs(s.x)), o.bound);
            if (~s.freezed)
                idx = find(1./(o.bound.*o.bound) < s.ham.n * s.ham.barrier.extraHessian);
                if ~isempty(idx)
                    s.ham.barrier.extraHessian = 0.1./(s.ham.n * o.bound.*o.bound);
                    s.ham.move(s.x, true); % Remember to update the cache used in barrier
                    s.v = s.ham.resample(s.x, zeros(size(s.x)));
                    s.log('DynamicRegularizer:set_bound', 'The bound of %i coordinates are changed.\n', length(idx));
                end
            end
        end
    end
end