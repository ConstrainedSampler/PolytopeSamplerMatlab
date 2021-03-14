function o = Solver(A, precision)
    if numel(A) == 0
        o = MatlabSolver(A);
        return;
    end
    
    o = MexSolver(A, precision);
end