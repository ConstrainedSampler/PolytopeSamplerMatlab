initSampler

N = 1000; % number of experiment

warning('off', 'uniformtest:size');
warning('off', 'stats:adtest:OutOfRangePLow');
warning('off', 'stats:adtest:OutOfRangePHigh');
warning('off', 'stats:adtest:SmallSampleSize');
warning('off', 'stats:adtest:NotEnoughData');

P = struct; d = 100;
P.Aineq = ones(1, d);
P.bineq = 1;
P.ub = ones(d, 1);
P.lb = zeros(d, 1);
p = zeros(N,1);
parfor i = 1:N
    opts = default_options();
    opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicWeight', 'DynamicStepSize'};
    opts.seed = i;
    o = sample(P, 1000, opts);
    p(i) = uniformtest(o);
end

% Look by eye to see if this is uniform
histogram(p, 20)
[~, z] = kstest(norminv(p));
assert(z > 0.05);

%% understand the worst case
if (0)
    opts = default_options();
    opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicWeight', 'DynamicStepSize'};
    opts.seed = find(p == min(p));
    o = sample(P, 1000, opts);
    uniformtest(o, struct('toPlot', true));
end