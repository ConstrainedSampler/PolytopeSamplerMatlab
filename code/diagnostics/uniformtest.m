function [pVal] = uniformtest(x, P, dim, opts)
%p_value = uniformtest(x, plan, opts)
%compute the p-value for the radial distribution of x
%
%Assumption: x follows uniform distribution in the interior of some convex set.
%
%Input:
% x - the samples object outputted by sample
% sampler - outputted by sample.
% opts - a structure for options with the following properties
%  toThin - extract independent samples from x (default: true)
%  toPlot - plot the radial distribution (to 1/(dim-1)) power
%
%Output:
% p_value - the p value of whether the empirical radial distribution follows the 
%           distribution r^(dim-1). We use Anderson-Darling and Kolmogorov-Smirnov tests here.

if ~exist('opts', 'var'), opts = struct; end
defaults.toThin = true;
defaults.toPlot = false;
defaults.tol = 1e-8;
opts = setfield(defaults, opts);
p = x(:,1);

if opts.toThin, x = thin_samples(x); end

if size(x,2) < 6
    warning('PolytopeSampler:uniformtest:size', 'Sample size must be at least 6.');
end

if isfield(P, 'df') && ~isempty(P.df)
    warning('PolytopeSampler:uniformtest:nonUniform', 'The density of the distribution should be uniform, namely, df = 0.');
end

K = size(x,2);

unif_vals = zeros(K,1);
for i=1:K
    this_x = x(:,i);
    u = this_x-p;
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

if opts.toPlot
    figure;
    cdfplot(unif_vals);
    hold on;
    plot(0:0.01:1, 0:0.01:1, '.')
end

try
    [~,pVal] = adtest(norminv(unif_vals));
catch
    pVal = 1;
end

end