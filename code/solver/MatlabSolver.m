classdef MatlabSolver < handle
    properties (GetAccess = public)
        % the constraint matrix
        A
        w = NaN
        w_solve = NaN
        
        % privates
        initialized = false
        L
    end
    
    methods
        % precision is either double or doubledouble
        function o = MatlabSolver(A)
            count = symbfact(A,'row');
            n = size(A,1);
            if sum(count.^2) > n^3 / 100
                
                A = sparse(A);
                H = A * A' + speye(n);
                t1 = timeit(@() chol(H, 'lower'));
                H = full(H);
                t2 = timeit(@() chol(H, 'lower'));
                
                if t1 > t2
                    A = full(A);
                end
            end
            
            o.A = A;
        end
        
        function setScale(o, w)
            o.initialized = true;
            
            if ~all(w == o.w) || (numel(w) ~= numel(o.w))
                H = o.A * (spdiag(w) * o.A');
                [o.L, p] = chol(H, 'lower');
                
                if p ~= 0
                    r = eps * (1 + sum(abs(H),1)');
                    while p ~= 0
                        [o.L, p] = chol(H + spdiag(r), 'lower');
                        r = r * 4 + eps;
                        if (any(r > 1e64))
                            error('Numerical error')
                        end
                    end
                end
            end
            
            o.w = w;
            o.w_solve = w;
        end
        
        function L = diagL(o)
            assert(o.initialized);
            
            L = diag(o.L);
        end
        
        function r = logdet(o)
            assert(o.initialized);
            
            r = 2 * sum(log(diag(o.L)));
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            assert(o.initialized);
            
            % TODO: Support non-JL
            % V = full((o.L\(o.A .* sqrt(o.w)'))');
            % sigma = 1-sum(V.^2,2);
            if nargin == 1 || nSketch == 0
                nSketch = 32;
            end
            
            JLdir = sign(rand(size(o.A,1), nSketch)-0.5);
            V = (o.A' * (o.L'\JLdir)) .* sqrt(o.w);
            sigma = 1 - sum(V.^2,2) / nSketch;
        end
        
        function y = approxSolve(o, b)
            assert(o.initialized);
            
            y = o.L'\(o.L\b);
        end
        
        function y = AwAt(o, b)
            assert(o.initialized);
            
            y = o.A * (o.w_solve .* (o.A' * b));
        end
        
        function x = solve(o, b, w, x0)
            assert(o.initialized);
            
            if nargin < 4, x0 = zeros(size(b)); end
            if nargin >= 3, o.w_solve = w; end
            x = batch_pcg(x0, b, o, 1e-12, 20);
        end
    end
end