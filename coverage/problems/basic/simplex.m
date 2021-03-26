function P = simplex(dim)
    P = struct;
    P.Aeq = ones(1,dim);
    P.beq = 1;
    P.lb = zeros(dim,1);
end