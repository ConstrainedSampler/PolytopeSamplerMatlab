classdef Hamiltonian < handle
    % H(x, v) = U(x) + K(x,v)
    % where
    % U(x) = f(x) + 1/2 (log det g + log det A g^-1 A')
    % K(x,v) = 1/2 v' (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v
    
    properties
        A       % constraint matrix A
        b       % constraint vector b
        f       % the objective function and its derivatives in the original space
        barrier % TwoSidedBarrier
        P       % domain
        lsc
        df
        ddf
        solver
        opts
    end
    
    % Dependent Varialbes
    properties
        m
        n
        k
        x
        fx % f(Tx + y)
        dfx % T' * f'(Tx + y)
        hess
        prepared
        last_dUdx = []
    end
    
    methods
        function o = Hamiltonian(P, opts)
            P.set_vdim(2);
            
            o.P = P;
            o.A = P.A;
            o.b = P.b;
            o.f = P.f;
            o.df = P.df;
            o.ddf = P.ddf;
            o.opts = opts;
            o.m = size(P.A,1);
            o.n = size(P.A,2);
            o.k = opts.simdLen;
            o.x = randn(o.k, o.n);
            o.barrier = P.barrier;
            o.solver = Solver(P.A, opts.solverThreshold, o.k);
            o.prepared = false;
        end
        
        % when we prepare 
        function prepare(o, x)
            o.move(x);
            if ~o.prepared
                o.solver.setScale(1./o.hess);
                o.last_dUdx = [];
            end
            o.prepared = true;
        end
        
        % Resample v = g^{1/2} * N(0, I_d)
        function v = resample(o, x, v, momentum)
            if nargin == 3, momentum = 0; end
            o.move(x);
            sqrtHess = sqrt(o.hess);
            v = sqrt(momentum) * v + sqrt(1-momentum) * (sqrtHess .* randn(o.k, o.n));
        end
        
        % Compute H(x,v)
        function E = H(o, x, v)
            o.prepare(x)
            K = 0.5 * sum(v .* o.DK(x, v),2);
            U = 0.5 * (o.solver.logdet() + sum(log(o.hess),2));
            U = U + o.fx;
            E = U + K;
        end
        
        function dUdx = DU(o, x)
            o.move(x);
            if ~o.prepared || isempty(o.last_dUdx)
                o.prepare(x);
                o.lsc = o.solver.leverageScoreComplement(o.opts.nSketch);
                o.last_dUdx = o.barrier.tensor(x) .* o.lsc ./ (2*o.hess) + o.dfx;
            end
            dUdx = o.last_dUdx;
        end
        
        % Project to Ax = b
        function x = project(o, x)
            o.move(x);
            % col vector: x = x + step * (o.A' * o.solver.approxSolve(o.b - o.A*x))./o.hess;
            x = x + (o.solver.approxSolve(o.b' - x*o.A') * o.A)./o.hess;
        end
        
        % Compute dK/dv = (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v and
        %         dK/dx = -Dg[dK/dv,dK/dv]/2
        function [dKdv, dKdx] = DK(o, x, v)
            o.move(x);
            invHessV = v./o.hess;
            % col vector: dKdv = invHessV - (o.A' * o.solver.solve(o.A * invHessV))./o.hess;
            dKdv = invHessV - (o.solver.solve(invHessV * o.A') * o.A)./o.hess;
            if nargout > 1
                dKdx = -o.barrier.quadratic_form_gradient(x, dKdv)/2;
            end
        end
        
        % Approximate dK/dv = (g^-1 - g^-1 A'(A g^-1 A')^-1 A g^-1) v and
        %             dK/dx = -Dg[dK/dv,dK/dv]/2
        function [dKdv, dKdx, nu] = approxDK(o, x, v, nu)
            o.move(x);
            % col vector: dUdv_b = o.A * ((v - o.A' * nu)./o.hess);
            dUdv_b = ((v - nu * o.A)./o.hess) * o.A';
            nu = nu + o.solver.approxSolve(dUdv_b);
            % col vector: dKdv = (v - o.A' * nu)./o.hess;
            dKdv = (v - nu * o.A)./o.hess;
            dKdx = -o.barrier.quadratic_form_gradient(x, dKdv)/2;
        end
        
        function t = step_size(o, x, dx)
            t1 = o.barrier.step_size(x, dx);
            t2 = 1 / max(sqrt(o.barrier.hessian_norm(x, dx)),[],2);
            t = min(t1,t2);
        end
        
        % Test if the values of x and v are valid and if x is feasible
        function r = feasible(o, x, v)
            r = ~any(isnan(x),2) & o.barrier.feasible(x);
            
            if nargin == 3
                r = r & ~any(isnan(v), 2);
            end
        end
        
        function r = v_norm(o, x, dv)
            o.move(x);
            r = sum((dv .* dv) ./ o.hess,2);
        end
        
        function r = x_norm(o, x, dx)
            o.move(x);
            r = sum((dx .* dx) .* o.hess,2);
        end
        
        function move(o, x, forceUpdate)
            if nargin == 2, forceUpdate = false; end
            if ~all(size(o.x) == size(x)), return; end
            if ~forceUpdate && all(o.x == x, 'all'), return; end
            
            
            o.x = x;
            [o.fx, o.dfx, ddfx] = o.P.f_oracle(x);
            o.hess = o.barrier.hessian(x) + ddfx;
            o.prepared = false;
        end
    end
end