classdef Hamiltonian < handle
    % H(x, v) = U(x) + K(x,v)
    % where
    % U(x) = f(x) + 1/2 (log det g + log det A g^-1 A')
    % K(x,v) = 1/2 v' (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v
    
    properties
        A   % constraint matrix A
        b   % constraint vector b
        T	% transformation T of the domain
        y	% shift of the domain
        f	% the objective function and its derivatives in the original space
        barrier % TwoSidedBarrier
        df
        ddf
        dddf
        solver
        crudeSolver
        accSolver
        precision
        accuracy
        opts
    end
    
    % Dependent Varialbes
    properties
        m
        n
        T2
        T3
        x
        hess
        prepared
        last_dUdx = []
        last_lsc = []
    end
    
    methods
        function o = Hamiltonian(P, opts)
            m = size(P.A,1); n = size(P.A,2);
            assert(all(size(P.b) == [m 1]));
            assert(P.barrier.n == n);
            
            o.A = P.A;
            o.b = P.b;
            o.f = P.f;
            o.df = P.df;
            o.ddf = P.ddf;
            o.dddf = P.dddf;
            o.opts = opts;
            o.m = m;
            o.n = n;
            o.x = randn(n, 1);
            o.barrier = P.barrier;
            o.crudeSolver = Solver(P.A, 'double');
            o.accSolver = Solver(P.A, 'doubledouble');
            o.solver = o.crudeSolver;
            o.prepared = false;
            o.precision = 'double';
            
            % assume each row of T contains at most 1 non-zero
            o.T = P.T;
            o.y = P.y;
            assert(max(full(sum(o.T~=0,2))) <= 1);
            o.T2 = o.T.^2;
            o.T3 = o.T.*o.T2;
            
            if ~isempty(o.df) && isfloat(o.df)
                o.df = o.T'*o.df;
            end
            if ~isempty(o.ddf) && isfloat(o.ddf)
                o.ddf = o.T2'*o.ddf;
            end
            if ~isempty(o.dddf) && isfloat(o.dddf)
                o.dddf = o.T3'*o.dddf;
            end
        end
        
        % when we prepare 
        function prepare(o, x)
            o.move(x);
            if ~o.prepared
                if o.opts.checkPrecision || strcmp(o.precision, 'doubledouble')
                    % first try if double precision is accurate enough
                    o.crudeSolver.setScale(1./o.hess);
                    z = 1 - o.crudeSolver.leverageScoreComplement(1);
                    o.accuracy = abs(sum(z)/(size(o.A,1)+eps) - 1); % + eps is to avoid 0/0 case
                    if o.accuracy < o.opts.crudeSolverThreshold
                        o.solver = o.crudeSolver;
                        o.precision = 'double';
                    else
                        o.accSolver.setScale(1./o.hess);
                        o.solver = o.accSolver;
                        o.precision = 'doubledouble';
                    end
                else
                    o.crudeSolver.setScale(1./o.hess);
                    o.accuracy = 0;
                end
                o.last_dUdx = [];
            end
            o.prepared = true;
        end
        
        % Resample v = g^{1/2} * N(0, I_d)
        function v = resample(o, x, v, momentum)
            if nargin == 3, momentum = 0; end
            o.move(x);
            sqrtHess = sqrt(o.hess);
            v = sqrt(momentum) * v + sqrt(1-momentum) * (sqrtHess .* randn(o.n,1));
        end
        
        % Compute H(x,v)
        function E = H(o, x, v)
            o.prepare(x)
            K = 0.5 * v' * o.DK(x, v);
            U = 0.5 * (o.solver.logdet() + sum(log(o.hess)));
            if ~isempty(o.df)
                if isfloat(o.df)
                    U = U + o.df' * x;
                else
                    U = U + o.f(o.T * x + o.y);
                end
            end
            
            E = U + K;
        end
        
        function dUdx = DU(o, x)
            o.move(x);
            if ~o.prepared || isempty(o.last_dUdx)
                o.prepare(x);
                o.last_lsc = o.solver.leverageScoreComplement(o.opts.nSketch);
                o.last_dUdx = o.barrier.tensor(x) .* o.last_lsc ./ (2*o.hess);
            end
            dUdx = o.last_dUdx;
            
            if ~isempty(o.df)
                if isfloat(o.df)
                    dUdx = dUdx + o.df;
                else
                    dUdx = dUdx + o.T' * o.df(o.T * x + o.y);
                end
            end
        end
        
        % Project to Ax = b
        function x = project(o, x, step)
            o.move(x);
            x = x + step * (o.A' * o.solver.approxSolve(o.b - o.A*x))./o.hess;
        end
        
        % Compute dK/dv = (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v and
        %         dK/dx = -Dg[dK/dv,dK/dv]/2
        function [dKdv, dKdx] = DK(o, x, v)
            o.move(x);
            invHessV = v./o.hess;
            dKdv = invHessV - (o.A' * o.solver.solve(o.A * invHessV))./o.hess;
            if nargout > 1
                dKdx = -o.barrier.quadratic_form_gradient(x, dKdv)/2;
            end
        end
        
        % Approximate dK/dv = (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v and
        %             dK/dx = -Dg[dK/dv,dK/dv]/2
        function [dKdv, dKdx, nu] = approxDK(o, x, v, nu)
            o.move(x);
            dUdv_b = o.A * ((v - o.A' * nu)./o.hess);
            nu = nu + o.solver.approxSolve(dUdv_b);
            dKdv = (v - o.A' * nu)./o.hess;
            dKdx = -o.barrier.quadratic_form_gradient(x, dKdv)/2;
        end
        
        function t = step_size(o, x, dx)
            t1 = o.barrier.step_size(x, dx);
            t2 = 1 / max(sqrt(o.barrier.hessian_norm(x, dx)));
            t = min(t1,t2);
        end
        
        % Test if the values of x and v are valid and if x is feasible
        function r = feasible(o, x, v)
            r = all(size(x) == [o.n 1]) && all(size(v) == [o.n 1]) && ...
                ~any(isnan(x)) && ~any(isnan(v)) && o.barrier.feasible(x);
        end
        
        function r = v_norm(o, x, dv)
            o.move(x);
            r = dv' * (dv./o.hess);
        end
        
        function r = x_norm(o, x, dx)
            o.move(x);
            r = dx' * (o.hess .* dx);
        end
        
        function move(o, x, forceUpdate)
            if nargin == 2, forceUpdate = false; end
            if ~forceUpdate && all(o.x == x), return; end
            
            o.x = x;
            h = o.barrier.hessian(x);
            
            if ~isempty(o.ddf)
                if isfloat(o.ddf)
                    h = h + o.ddf;
                else
                    h = h + o.T2' * o.ddf(o.T * x + o.y);
                end
            end
            
            o.hess = h;
            o.prepared = false;
        end
    end
end