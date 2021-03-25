classdef MultiMatlabSolver2 < handle
    properties (GetAccess = public)
        solvers
        k
        m
        n
        w = NaN
        accuracy
        uid
        solver
        w_solve
        A
    end
    
    methods
        % precision is either double or doubledouble
        function o = MultiMatlabSolver2(A, precision, k)
            o.k = k;
            o.m = size(A,2);
            o.n = size(A,1);
            for i = 1:o.k
                o.solvers{i} = MexSolver(A, precision, 0);
            end
            
            o.solver = str2func(['PackedChol' num2str(k)]);
            o.uid = o.solver('init', uint64(randi(2^32-1,'uint32')), A);
            o.solver('setAccuracyTarget', o.uid, precision);
            o.A = A;
        end
        
        function delete(o)
            o.solver('delete', o.uid);
        end
        
        function setScale(o, w)
            o.accuracy = zeros(o.k, 1);
            for i = 1:o.k
                o.solvers{i}.setScale(reshape(w(i,:), [o.m, 1]));
                o.accuracy(i) = o.solvers{i}.accuracy;
            end
            
            
            if ~all(w == o.w, 'all')
                o.accuracy = o.solver('decompose', o.uid, w);
                o.w = w;
            end
            
            o.w_solve = w;
        end
        
        function L = diagL(o)
            L = zeros(o.k, o.n);
            for i = 1:o.k
                L(i, :) = o.solvers{i}.diagL();
            end
            
            L2 = o.solver('diagL', o.uid);
            assert(~(norm(L-L2, 'fro') > 1e-4 * norm(L, 'fro')));
        end
        
        function r = logdet(o)
            r = zeros(o.k, 1);
            for i = 1:o.k
                r(i) = o.solvers{i}.logdet();
            end
            
            r2 = o.solver('logdet', o.uid);
            assert(~(norm(r-r2, 'fro') > 1e-4 * norm(r, 'fro')));
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            if nargin == 1, nSketch = 0; end
            
            sigma = zeros(o.k, o.m);
            for i = 1:o.k
                sigma(i,:) = o.solvers{i}.leverageScoreComplement(nSketch);
            end
            
            sigma2 = o.solver('leverageScoreComplement', o.uid, nSketch);
            assert(~(norm(sigma-sigma2, 'fro') > 1e-4 * norm(sigma, 'fro')));
        end
        
        function y = approxSolve(o, b)
            y = zeros(o.k, o.n);
            for i = 1:o.k
                y(i,:) = o.solvers{i}.approxSolve(b(i,:)');
            end
            
            y2 = o.solver('solve', o.uid, b);
            assert(~(norm(y-y2, 'fro') > 1e-4 * norm(y, 'fro')));
        end
        
        function x = solve(o, b, w, x0)
            if nargin < 4
                x0 = zeros([o.k, o.n]);
            end
            
            if nargin < 3, w = o.w; end
            
            x = zeros(o.k, o.n);
            for i = 1:o.k
                x(i,:,:) = o.solvers{i}.solve(b(i,:)', w(i,:)', x0(i,:)');
            end
            
            x2 = batch_pcg(x0, b, o, 1e-12, 20);
            assert(~(norm(x-x2, 'fro') > 1e-4 * norm(x, 'fro')));
        end
        
        function counts = getDecomposeCount(o)
            counts = zeros(o.k+1, 1);
            for i = 1:o.k
                c = o.solvers{i}.getDecomposeCount();
                counts(i) = c(1);
                counts(end) = c(2);
            end
        end
        
        % used only in batch_pcg
        function y = AwAt(o, b)
            if (o.k == 0)
                y = o.A * (o.w_solve .* (o.A' * b));
            else
                y = ((b * o.A) .* o.w_solve) * o.A';
            end
        end
    end
end