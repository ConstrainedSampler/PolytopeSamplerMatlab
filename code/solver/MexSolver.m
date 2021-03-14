classdef MexSolver < handle
    properties
        % the constraint matrix
        A
        w = NaN
        w_solve = NaN
        
        % private
        uid
        uid2
        initialized = false
        precision
    end
    
    methods
        % precision is either double or doubledouble
        function o = MexSolver(A, precision)
            if nargin < 2, precision = 'double'; end
            
            o.A = A;
            o.uid = PackedChol('init', uint64(randi(2^32-1,'uint32')), A);
            if (strcmp(precision,'double'))
                PackedChol('setAccuracyTarget', o.uid, 1e-6);
            elseif (strcmp(precision,'doubledouble'))
                PackedChol('setAccuracyTarget', o.uid, 0.0);
            else
                error('Unsupported precision mode');
            end
        end
        
        function delete(o)
            PackedChol('delete', o.uid);
        end
        
        function setScale(o, w)
            o.initialized = true;
            if ~all(w == o.w) || (numel(w) ~= numel(o.w))
                PackedChol('decompose', o.uid, w);
            end
            
            o.w = w;
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
        
        function y = AwAt(o, b)
            assert(o.initialized);
            
            y = o.A * (o.w_solve .* (o.A' * b));
        end
    end
end