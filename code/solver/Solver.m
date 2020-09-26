function o = Solver(A, precision)
    if numel(A) == 0
        o = MatlabSolver(A);
        return;
    end
    
    if nargin == 1 || strcmp(precision, 'double')
        o1 = MexSolver(A, 'double');
        o2 = MatlabSolver(A);
        m = size(o1.A,2);
        
        t1 = timeit(@() o1.setScale(rand(m,1) + 1));
        t2 = timeit(@() o2.setScale(rand(m,1) + 1));
        
        
        if t1 < t2
            o = o1;
        else
            o = o2;
        end
    else
        o = MexSolver(A, precision);
    end
end