classdef TwoSidedBarrier < handle
% The barrier for the domain {lu <= x <= ub}.
% phi(x) = - sum log(x - lb) - sum log(ub - x)
    properties (SetAccess = private)
        upperIdx % upper only index 
        lowerIdx % lower only index
        freeIdx % free index
        center  % a feasible point
        n % number of variables
        ub % ub
        lb % lb
    end
    
    properties
        extraHessian % the extra constant Hessian added for free constraints
    end
    
    methods
        function o = TwoSidedBarrier(lb, ub)
            o.update(lb, ub);
        end
        
        function o = update(o, lb, ub)
            o.n = length(lb);
            assert(all(size(lb) == [o.n 1]));
            assert(all(size(ub) == [o.n 1]));
            assert(all(lb < ub));
            
            o.lb = lb;
            o.ub = ub;
            o.upperIdx = find(o.lb == -Inf);
            o.lowerIdx = find(o.ub == Inf);
            o.freeIdx = find((o.lb == -Inf) & (o.ub == Inf));
            
            c = (o.ub+o.lb)/2;
            c(o.lowerIdx) = o.lb(o.lowerIdx) + 1e6;
            c(o.upperIdx) = o.ub(o.upperIdx) - 1e6;
            c(o.freeIdx) = 0;
            o.center = c;
        end
        
        function r = feasible(o, x)
            % r = feasible(o, x)
            % output if x is feasible
            
            r = all((x > o.lb) & (x < o.ub));
        end
        
        function t = step_size(o, x, v)
            % t = stepsize(o, x, v)
            % output the maximum step size from x with direction v
            
            % check positive direction
            posIdx = v > 0;
            t = min([+1e40;(o.ub(posIdx) - x(posIdx))./v(posIdx)]); %1e40 is max for single
            
            % check negative direction
            negIdx = v < 0;
            t = min([t;(o.lb(negIdx) - x(negIdx))./v(negIdx)]);
        end
        
        function [A, b] = boundary(o, x)
            % [A, b] = boundary(o, x)
            % output the normal at the boundary around x for each barrier
            
            c = o.center;
            
            b = o.ub;
            b(x<c) = -o.lb(x<c);
            
            A = ones(size(x));
            A(x<c) = -A(x<c);
            
            A = spdiag(A);
        end
        
        function grad = gradient(o, x)
            grad = 1./(o.ub-x) - 1./(x-o.lb);
        end
        
        function d = hessian(o, x)
            d = 1./((x-o.lb).*(x-o.lb)) + 1./((o.ub-x).*(o.ub-x)) + o.extraHessian;
        end
        
        function t = tensor(o, x)
            t = -2*(1./((x-o.lb).*(x-o.lb).*(x-o.lb)) - 1./((o.ub-x).*(o.ub-x).*(o.ub-x)));
        end
        
        function v = quadratic_form_gradient(o, x, u)
            % v = quadratic_form_gradient(o, x, u)
            % output the dense vector - grad of sum_i u_i^T (hess phi(x)) u_i
            % where each col of u is one vector
            
            t = o.tensor(x);
            v = sum(u.^2,2).*t;
        end
        
        function v = logdet_gradient(o, x) 
            % output the dense vector - gradient of log det (hess phi(x))
            
            d = o.hessian(x);
            t = o.tensor(x);
            v = t./d;
        end
        
        % TODO: Add explanation
        function v = boundary_distance(o, x)
            v = abs(min(x-o.lb, o.ub-x));
        end
        
        % TODO: Add explanation
        function u = hessian_norm(o, x, v)
            d = o.hessian(x);
            u = (v .* v) .* d;
        end
    end
end