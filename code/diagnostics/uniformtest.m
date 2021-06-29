function [pVal] = uniformtest(o, opts)
%p_value = uniformtest(o, opts)
%Compute the p-value for the radial distribution of o.samples
%
%Assumption: o.samples follows uniform distribution in the interior of some convex set.
%
%Input:
% o - The samples object outputted by sample
% opts - A structure for options with the following properties
%  toPlot - Plot the radial distribution (to 1/(dim-1)) power (default = false)
%  thinningFactor - The factor we further sample the o.samples (default = 5)
%  tol - The amount of infeasibility we allow for each constraints (default = 1e-8)
%  
%
%Output:
% pVal - The p value of whether the empirical radial distribution follows the
%        distribution r^(dim-1).
%        We use Anderson-Darling and Kolmogorov-Smirnov tests here.

if nargin == 1, opts = struct; end

if isempty(ver('stats'))
    warning('uniformtest requires the Statistics and Machine Learning Toolbox');
    return
end
assert(~o.opts.rawOutput, 'This function does not support raw output');

defaults.toPlot = false;
defaults.thinningFactor = 5;
defaults.tol = 1e-8;
opts = setfield(defaults, opts);

% Since adtest and kstest are for independent samples,
% we further sample the samples to make sure they are independent enough.
x = o.samples(:, 1:opts.thinningFactor:end);

% The Hausdorff dimension of the domain
dim = o.problem.n - size(o.problem.A, 1);

if size(x,2) < 20
    warning('uniformtest:size', 'Effective sample size should be at least %i.', opts.thinningFactor * 20);
    if size(x,2) < 5
        pVal = 0;
        return;
    end
end

p = x(:,1);
x = x(:,2:end);

P = o.problem.original;
if ~isempty(P.df)
    warning('PolytopeSampler:uniformtest:nonUniform', 'This test only works for uniform distributions, namely, df = 0.');
end

N = size(x,2);
unif_vals = zeros(N,1);
for i = 1:N
    x_i = x(:,i);
    u = x_i - p;
    
    % compute r = ||x_i||_K
    
    % check the constraint x <= P.ub
    posIdx = u > opts.tol;
    r1 = min((P.ub(posIdx) - x_i(posIdx))./(u(posIdx)));
    
    % check the constraint x >= P.lb
    negIdx = u < -opts.tol;
    r2 = min((P.lb(negIdx) - x_i(negIdx))./(u(negIdx)));
    
    % check the constraint P.Aineq <= P.bineq
    Ax = P.Aineq * x_i; Au = P.Aineq * u;
    posAIdx = (Au > opts.tol);
    r3 = min((P.bineq(posAIdx) - Ax(posAIdx))./(Au(posAIdx)));
    
    r = min([r1, r2, r3]);
    unif_vals(i) = (1 / (1+r))^dim;
end

% To ensure pval1 and pval2 are independent,
% we run adtest and kstest on disjoint subset
try
    [~,pVal1] = adtest(norminv(unif_vals(1:floor(N/2))));
catch
    pVal1 = 0;
end

try
    [~,pVal2] = kstest(norminv(unif_vals((floor(N/2)+1):end)));
catch
    pVal2 = 0;
end

% We merge two indepdent p value into one p value
z = pVal1 * pVal2;
pVal = z - z * log(z);

if opts.toPlot
    figure;
    cdfplot(unif_vals);
    hold on;
    plot(0:0.01:1, 0:0.01:1, '.')
    title(sprintf('Empirical CDF of the radial distribution (pval = %.4f)', pVal))
end