classdef TwoSidedBarrier < handle
% The barrier for the domain {lu <= x <= ub}.
% phi(x) = - sum log(x - lb) - sum log(ub - x)
% if x is a matrix, each x is along the "dim"-th dim
    properties (SetAccess = private)
        upperIdx % upper only index 
        lowerIdx % lower only index
        freeIdx % free index
        center  % a feasible point
        n % number of variables
        ub % ub
        lb % lb
        opDim % each vector is along the "dim"-th dim
    end
    
    properties
        extraHessian % the extra constant Hessian added for free constraints
    end
    
    methods
        function o = TwoSidedBarrier(lb, ub, opDim)
            o.update(lb, ub);
            if nargin < 3, opDim = 1; end
            o.opDim = opDim;
        end
        
        function o = update(o, lb, ub)
            o.n = length(lb);
            assert(numel(lb) == o.n);
            assert(numel(ub) == o.n);
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
            r = all((x > o.lb) & (x < o.ub), o.opDim);
        end
        
        function t = step_size(o, x, v)
            % t = stepsize(o, x, v)
            % output the maximum step size from x with direction v
            
            max_step = 1e16; % largest step size
            
            % check positive direction
            posIdx = v > 0;
            t1 = min((o.ub(posIdx) - x(posIdx))./v(posIdx), [], o.opDim);
            if isempty(t1), t1 = max_step; end
            
            % check negative direction
            negIdx = v < 0;
            t2 = min((o.lb(negIdx) - x(negIdx))./v(negIdx), [], o.opDim);
            if isempty(t2), t2 = max_step; end
            
            t = min(min(t1, t2), max_step);
        end
        
        function [A, b] = boundary(o, x)
            % [A, b] = boundary(o, x)
            % output the normal at the boundary around x for each barrier
            % Assume: only 1 vector is given
            
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
            % output the -grad of u' (hess phi(x)) u
            
            t = o.tensor(x);
            v = u.*u.*t;
        end
        
        function v = logdet_gradient(o, x) 
            % output the gradient of log det (hess phi(x))
            % which is (hess phi(x))^-1 (grad hess phi (x))
            
            d = o.hessian(x);
            t = o.tensor(x);
            v = t./d;
        end
        
        function v = boundary_distance(o, x)
            % compute the distance of x and the closest boundary for each
            % constraint
            v = abs(min(x-o.lb, o.ub-x));
        end
        
        function u = hessian_norm(o, x, d)
            % compute d_i Hessian(x)_ii d_i for each i
            h = o.hessian(x);
            u = (d .* d) .* h;
        end
        
        function selectOpDim(o, opDim)
            o.opDim = opDim;
            if opDim == 1
                o.lb = reshape(o.lb, [o.n, 1]);
                o.ub = reshape(o.ub, [o.n, 1]);
                o.c = reshape(o.c, [o.n, 1]);
            else
                o.lb = reshape(o.lb, [1, o.n]);
                o.ub = reshape(o.ub, [1, o.n]);
                o.center = reshape(o.center, [1, o.n]);
            end
        end
    end
end