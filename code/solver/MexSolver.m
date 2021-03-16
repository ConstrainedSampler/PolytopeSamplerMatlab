classdef MexSolver < handle
    properties
        % the constraint matrix
        A
        w = NaN
        w_solve = NaN
        
        % private
        uid
        initialized = false
        precision
    end
    
    % Okay, the simd version simply do it one by one
    % w is a matrix k by n
    % b is a matrix k by n by ??
    % b, w, x0
    
    methods
        % precision is either double or doubledouble
        function o = MexSolver(A, precision)
            o.A = A;
            o.uid = PackedChol('init', uint64(randi(2^32-1,'uint32')), A);
            PackedChol('setAccuracyTarget', o.uid, precision);
        end
        
        function delete(o)
            PackedChol('delete', o.uid);
        end
        
        function setScale(o, w)
            o.initialized = true;
            if ~all(w == o.w)
                PackedChol('decompose', o.uid, w);
                o.w = w;
            end
            
            o.w_solve = w;
        end
        
        function r = diagL(o)
            assert(o.initialized);
            
            r = PackedChol('diagL', o.uid);
        end
        
        function r = logdet(o)
            assert(o.initialized);
            
            r = PackedChol('logdet', o.uid);
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            assert(o.initialized);
            
            if nargin == 1, nSketch = 0; end
            sigma = PackedChol('leverageScoreComplement', o.uid, nSketch);
        end
        
        function counts = getDecomposeCount(o)
            counts = PackedChol('getDecomposeCount', o.uid);
        end
        
        function y2 = approxSolve(o, b)
            assert(o.initialized);
            
            y2 = PackedChol('solve', o.uid, b')';
        end
        
        function x = solve(o, b, w, x0)
            assert(o.initialized);
            
            if nargin < 4, x0 = zeros(size(b)); end
            if nargin >= 3, o.w_solve = w; end
            x = batch_pcg(x0, b, o, 1e-12, 20);
        end
        
        % used only in batch_pcg
        function y = AwAt(o, b)
            assert(o.initialized);
            
            y = o.A * (o.w_solve .* (o.A' * b));
        end
    end
end