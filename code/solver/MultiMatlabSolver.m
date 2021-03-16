classdef MultiMatlabSolver < handle
    properties (GetAccess = public)
        solvers
        k
        m
        n
        w
    end
    
    methods
        % precision is either double or doubledouble
        function o = MultiMatlabSolver(A, precision, k)
            o.k = k;
            o.m = size(A,2);
            o.n = size(A,1);
            for i = 1:o.k
                o.solvers{i} = MatlabSolver(A, precision);
            end
        end
        
        function setScale(o, w)
            w = reshape(w, [o.k, o.m]);
            o.w = w;
            for i = 1:o.k
                o.solvers{i}.setScale(reshape(w(i,:), [o.m, 1]));
            end
        end
        
        function L = diagL(o)
            L = zeros(o.k, o.n);
            for i = 1:o.k
                L(i, :) = o.solvers{i}.diagL();
            end
        end
        
        function r = logdet(o)
            r = zeros(o.k, 1);
            for i = 1:o.k
                r(i) = o.solvers{i}.logdet();
            end
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            if nargin == 1, nSketch = 0; end
            
            sigma = zeros(o.k, o.m);
            for i = 1:o.k
                sigma(i,:) = o.solvers{i}.leverageScoreComplement(nSketch);
            end
        end
        
        function y = approxSolve(o, b)
            l = numel(b)/(o.k*o.n);
            b = reshape(b, [o.k, o.n, l]);
            
            y = zeros(o.k, o.n, l);
            for i = 1:o.k
                b_ = reshape(b(i,:,:), [o.n, l]);
                y(i,:,:) = o.solvers{i}.approxSolve(b_);
            end
        end
        
        function x = solve(o, b, w, x0)
            l = numel(b)/(o.k*o.n);
            b = reshape(b, [o.k, o.n, l]);
            if nargin < 4
                x0 = zeros([o.k, o.n, l]);
            else
                x0 = reshape(x0, [o.k, o.n, l]);
            end
            if nargin < 3, w = o.w; end
            w = reshape(w, [o.k, o.m]);
            
            x = zeros(o.k, o.n, l);
            for i = 1:o.k
                b_ = reshape(b(i,:,:), [o.n, l]);
                x0_ = reshape(x0(i,:,:), [o.n, l]);
                w_ = reshape(w(i,:), [o.m, 1]);
                x(i,:,:) = o.solvers{i}.solve(b_, w_, x0_);
            end
        end
        
        function counts = getDecomposeCount(o)
            counts = zeros(o.k+1);
            for i = 1:o.k
                c = o.solvers{i}.getDecomposeCount();
                counts(i) = c(1);
                counts(end) = c(2);
            end
        end
    end
end