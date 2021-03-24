initSampler

for k = 2:5
    P = struct; d = 10^k;
    P.ub = ones(d, 1);
    P.lb = zeros(d, 1);

    opts = default_options();
    opts.simdLen = 1;
    opts.odeMethod = @gauss_legendre;
    opts.MixingTimeEstimator.startIter = 1000;
    o = sample(P, 10, opts);
end