classdef MatlabSolver < handle
    properties (GetAccess = public)
        % the constraint matrix
        A
        w = NaN
        w_solve = NaN
        
        % privates
        initialized = false
        L
        precision
        usingExact = false
        numExact
        accuracy
        exactSolver
    end
    
    methods
        % precision is either double or doubledouble
        function o = MatlabSolver(A, precision)
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
            o.precision = precision;
            o.numExact = zeros(2,1);
            o.exactSolver = MexSolver(A, 0.0, 0);
        end
        
        function err = estimateAccuracy(o)
            k = 4;
            JLdir = (rand(size(o.A,1), 2 * k)-0.5)*sqrt(12);
            V = (o.A' * (o.L'\JLdir)) .* sqrt(o.w);
            err = 0;
            for i = 1:k
                err = err + (sum(V(:,2*i-1).*V(:,2*i)) - sum(JLdir(:,2*i-1).*JLdir(:,2*i)))^2;
            end
            err = sqrt(err / k);
        end
        
        function setScale(o, w)
            o.initialized = true;
            
            if ~all(w == o.w)
                o.w = w;
                o.numExact(2) = o.numExact(2) + 1;
                o.usingExact = false;
                if (o.precision > 0.0)
                    H = o.A * (spdiag(w) * o.A');
                    [o.L, p] = chol(H, 'lower');
                    if p~= 0
                        o.usingExact = true;
                    else
                        o.accuracy = o.estimateAccuracy();
                        if (isnan(o.accuracy) || o.accuracy > o.precision)
                            o.usingExact = true;
                        end
                    end
                else
                    o.usingExact = true;
                end
                
                if (o.usingExact)
                    o.numExact(1) = o.numExact(1) + 1;
                    o.exactSolver.setScale(w);
                end
            end
            
            o.w_solve = w;
        end
        
        function L = diagL(o)
            assert(o.initialized);
            
            if (o.usingExact)
                L = o.exactSolver.diagL();
            else
                L = diag(o.L);
            end
        end
        
        function r = logdet(o)
            assert(o.initialized);
            
            if (o.usingExact)
                r = o.exactSolver.logdet();
            else
                r = 2 * sum(log(diag(o.L)));
            end
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            assert(o.initialized);
            
            if (o.usingExact)
                if nargin == 1
                    sigma = o.exactSolver.leverageScoreComplement();
                else
                    sigma = o.exactSolver.leverageScoreComplement(nSketch);
                end
            else
                % V = full((o.L\(o.A .* sqrt(o.w)'))');
                % sigma = 1-sum(V.^2,2);
                if nargin == 1 || nSketch == 0
                    nSketch = 32;
                end

                JLdir = sign(rand(size(o.A,1), nSketch)-0.5);
                V = (o.A' * (o.L'\JLdir)) .* sqrt(o.w);
                sigma = 1 - sum(V.^2,2) / nSketch;
            end
        end
        
        function y = approxSolve(o, b)
            assert(o.initialized);
            
            if (o.usingExact)
                y = o.exactSolver.approxSolve(b);
            else
                y = o.L'\(o.L\b);
            end
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
        
        function counts = getDecomposeCount(o)
            counts = o.numExact;
        end
    end
end