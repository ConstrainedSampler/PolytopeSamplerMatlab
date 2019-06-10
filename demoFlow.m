%% Demo to sample from a flow polytope
initSampler

%% Form the problem P
n = 1000; e = ones(n,1);
P = Problem;
% the variables are [x, y] where x in R^n, y in R^{n-1}

% y_i = x_i - x_{i-1}
P.Aeq = [spdiags([e -e], 0:1, n-1, n) ...
    spdiags(-e, 0, n-1, n-1)];
P.beq = zeros(n-1,1);

% x_1 = 0
P.Aeq(end+1,1) = 1;
P.beq(end+1) = 0;

% x_n = 0
P.Aeq(end+1,n) = 1;
P.beq(end+1) = 0;

% -10 sqrt(n) <= x <= 10 * sqrt(n), -1 <= y <= 1
P.lb = -ones(2*n-1,1);
P.ub = ones(2*n-1,1);
P.lb(1:n) = -10*sqrt(n);
P.ub(1:n) = 10*sqrt(n);

%% Samples from P
iter = 10000;

tic;
plan = prepare(P, struct('display', 1, 'recordInterval', 10));
out = sample(plan, iter);
t = toc;

%% Output the result
fprintf('Total time = %f sec\n', t)

[ess] = effectiveSampleSize(out.samples);
fprintf('Mixing time = %f iter\n', size(out.samples,2) / min(ess))

p = unifScaleTest(out, plan, struct('toPlot',1));
fprintf('p value = %f\n', p)

figure
plot(out.samples(1:n,end));
title('Random sample from a flow polytope');