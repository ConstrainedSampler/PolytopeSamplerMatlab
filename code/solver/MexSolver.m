classdef MexSolver < handle
    properties
        % the constraint matrix
        A
        w = NaN
        w_solve = NaN
        precision
		
        % private
        uid
        initialized = false
        accuracy
        solver
        k
    end
    
    methods (Static)
        function o = loadobj(s)
            s.uid = s.solver('init', uint64(randi(2^32-1,'uint32')), s.A);
			s.solver('setAccuracyTarget', s.uid, s.precision);
            if ~any(isnan(s.w))
        		w = s.w; s.w = NaN;
            	s.setScale(w);
            end
            o = s;
        end
    end
	
    methods
        % precision is either double or doubledouble
        function o = MexSolver(A, precision, k)
            o.A = A;
            o.k = k;
            o.solver = str2func(['PackedChol' num2str(k)]);
            o.uid = o.solver('init', uint64(randi(2^32-1,'uint32')), A);
            o.solver('setAccuracyTarget', o.uid, precision);
			o.precision = precision;
        end
        
		function b = saveobj(a)
			b = a;
			b.uid = [];
        end
        
        function delete(o)
			if ~isempty(o.uid)
				o.solver('delete', o.uid);
			end
        end
        
        function setScale(o, w)
            o.initialized = true;
            if ~all(w == o.w, 'all')
                o.accuracy = o.solver('decompose', o.uid, w);
                o.w = w;
            end
            
            o.w_solve = w;
        end
        
        function r = diagL(o)
            assert(o.initialized);
            
            r = o.solver('diagL', o.uid);
        end
        
        function r = logdet(o)
            assert(o.initialized);
            
            r = o.solver('logdet', o.uid);
        end
        
        % Note that sigma can be negative. 
        % We cannot truncate it to >=0 because we need unbias
        function sigma = leverageScoreComplement(o, nSketch)
            assert(o.initialized);
            
            if nargin == 1, nSketch = 0; end
            sigma = o.solver('leverageScoreComplement', o.uid, nSketch);
        end
        
        function counts = getDecomposeCount(o)
            counts = o.solver('getDecomposeCount', o.uid);
        end
        
        function y2 = approxSolve(o, b)
            assert(o.initialized);
            
            y2 = o.solver('solve', o.uid, b);
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
            if (o.k == 0)
                y = o.A * (o.w_solve .* (o.A' * b));
            else
                y = ((b * o.A) .* o.w_solve) * o.A';
            end
        end
    end
end