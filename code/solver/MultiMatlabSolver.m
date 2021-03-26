classdef MultiMatlabSolver < handle
    properties (GetAccess = public)
        solvers
        k
        m
        n
        w
        accuracy
    end
    
    methods
        function o = MultiMatlabSolver(A, precision, k)
            o.k = k;
            o.m = size(A,2);
            o.n = size(A,1);
            for i = 1:o.k
                o.solvers{i} = MexSolver(A, precision, 0);
            end
        end
        
        function setScale(o, w)
            o.w = w;
            o.accuracy = zeros(o.k, 1);
            for i = 1:o.k
                o.solvers{i}.setScale(reshape(w(i,:), [o.m, 1]));
                o.accuracy(i) = o.solvers{i}.accuracy;
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
            y = zeros(o.k, o.n);
            for i = 1:o.k
                y(i,:) = o.solvers{i}.approxSolve(b(i,:)');
            end
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
        end
        
        function counts = getDecomposeCount(o)
            counts = zeros(o.k+1, 1);
            for i = 1:o.k
                c = o.solvers{i}.getDecomposeCount();
                counts(i) = c(1);
                counts(end) = c(2);
            end
        end
    end
end