function P = box(dim)
    P = struct;
    P.lb = -0.5*ones(dim,1);
    P.ub = 0.5*ones(dim,1);
end