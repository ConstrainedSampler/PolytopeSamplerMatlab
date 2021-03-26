function o = Solver(A, precision, k)
    if nargin < 2, precision = 'double'; end
    if nargin < 3, k = 0; end
    
    if (strcmp(precision,'double'))
        precision = 1e-6;
    elseif (strcmp(precision,'doubledouble'))
        precision = 0.0;
    elseif ~isfloat(precision)
        error('Unsupported precision mode');
    end

    o = MexSolver(A, precision, k);
    %if nargin < 3
    %    o = MatlabSolver(A, precision);
    %else
    %    o = MultiMatlabSolver(A, precision, k);
    %end
end