function [pVal] = uniformtest(o, opts)
%p_value = uniformtest(o, opts)
%compute the p-value for the radial distribution of o.samples
%
%Assumption: o.samples follows uniform distribution in the interior of some convex set.
%
%Input:
% o - the samples object outputted by sample
% opts - a structure for options with the following properties
%  toPlot - plot the radial distribution (to 1/(dim-1)) power
%
%Output:
% pVal - the p value of whether the empirical radial distribution follows the 
%        distribution r^(dim-1). We use Anderson-Darling and Kolmogorov-Smirnov tests here.

if ~exist('opts', 'var'), opts = struct; end

if ~exist('adtest') || ~exist('kstest')
    error('This function requires the Statistics and Machine Learning Toolbox');
end

defaults.toPlot = false;
defaults.tol = 1e-8;
opts = setfield(defaults, opts);
P = o.polytope.originalProblem;
x = thin_samples(o.samples);
s = size(x); l = prod(s(1:end-1));
x = reshape(x, [l, s(end)]);
dim = o.polytope.n - size(o.polytope.A, 1);

if size(x,2) < 10
    warning('uniformtest:size', 'Effective sample size should be at least 10.');
end

p = x(:,1);
x = x(:,2:end);

if ~isempty(P.df)
    warning('PolytopeSampler:uniformtest:nonUniform', 'The density of the distribution should be uniform, namely, df = 0.');
end

K = size(x,2);

unif_vals = zeros(K,1);
for i=1:K
    this_x = x(:,i);
    u = this_x - p;
    % P.Aineq <= P.bineq
    % P.lb <= x <= P.ub
    
    posIdx = u > opts.tol; % x <= P.ub
    r = min([+1e40;(P.ub(posIdx) - this_x(posIdx))./(u(posIdx))]);
    negIdx = u < -opts.tol; % P.lb <= x
    r = min([r;(P.lb(negIdx) - this_x(negIdx))./(u(negIdx))]);
    
    Ax = P.Aineq * this_x;
    Au = P.Aineq * u;
    posAIdx = (Au > 1e-8);
    r = min([r;(P.bineq(posAIdx) - Ax(posAIdx))./(Au(posAIdx))]);
    
    unif_vals(i) = (1 / (1+r))^dim;
end

try
    [~,pVal1] = adtest(norminv(unif_vals));
catch
    pVal1 = 1;
end

try
    [~,pVal2] = kstest(norminv(unif_vals));
catch
    pVal2 = 1;
end

z = pVal1 * pVal2;
pVal = z - z * log(z);

if opts.toPlot
    figure;
    cdfplot(unif_vals);
    hold on;
    plot(0:0.01:1, 0:0.01:1, '.')
    title(sprintf('Empirical CDF of the radial distribution (pval = %.4f)', pVal))
end

end