function [pVal] = distribution_test(o, opts)
%p_value = distribution_test(o, opts)
%Compute the p-value for the radial distribution of o.samples
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
opts = Setfield(defaults, opts);

% Since adtest and kstest are for independent samples,
% we further sample the samples to make sure they are independent enough.
x = o.samples(:, 1:opts.thinningFactor:end);

% The Hausdorff dimension of the domain
dim = o.problem.n - size(o.problem.A, 1);

if size(x,2) < 50
   warning('uniformtest:size', 'Effective sample size should be at least %i.', opts.thinningFactor * 50);
   if size(x,2) < 10
      pVal = 0;
      return;
   end
end

p = x(:,1);
x = x(:,2:end);

P = o.problem.original;
if isempty(P.df)
   P.df = zeros(size(x,1),1);
end
vectorMode = isfloat(P.df);

if (~vectorMode)
   pt = P.f(p) + 1;
   dim = dim + numel(pt);
end

N = size(x,2);
unif_vals = zeros(N,1);

pLbS = P.lb - p;
pUbS = P.ub - p;
pIneqS = P.bineq - P.Aineq * p;
for i = 1:N
   x_i = x(:,i);
   u = x_i - p;
   
   % compute r be largest scalar st p + r u is feasible
   
   % check the constraint x <= P.ub
   posIdx = u > opts.tol;
   r1 = min(pUbS(posIdx)./u(posIdx));
   
   % check the constraint x >= P.lb
   negIdx = u < -opts.tol;
   r2 = min(pLbS(negIdx)./(u(negIdx)));
   
   % check the constraint P.Aineq * x <= P.bineq
   Au = P.Aineq * u;
   posAIdx = (Au > opts.tol);
   r3 = min(pIneqS(posAIdx)./(Au(posAIdx)));
   
   r = min([r1, r2, r3]);
   if (vectorMode)
      a = u' * P.df;
   else
      xt = P.f(x_i);
      xt = xt + exprnd(1,size(xt,1),size(xt,2));
      ut = xt - pt;
      g = @(t) P.f(p + t * u) - (pt + t * ut);
      r4 = binary_search(g, 1, r, 1e-12);
      r = min(r, r4);
      
      a = sum(ut);
   end
   
   if (a * r > dim)
      v = gammainc(a, dim) / gammainc(a * r, dim);
   else
      v = scaledgammainc(a, dim) / scaledgammainc(a * r, dim);
      v = v * exp(a * (r-1) - dim * log(r));
   end
   unif_vals(i) = real(v);
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
end

% find maximum t such that g(x) <= 0
% assume g(a) <= 0
function t = binary_search(g, a, b, tol)
   assert(all(g(a) <= 0) && b >= a);
   if all(g(b)<=0)
      t = b;
   else
      while (b > a + tol)
         m = (b+a)/2;
         if all(g(m)<=0)
            a = m;
         else
            b = m;
         end
      end
      t = (a+b)/2;
   end
end

function t = scaledgammainc(a, d)
assert(a <= d);
if (a >= -d-20)
   t = gammainc(a, d, 'scaledlower');
else
   t = (d)/(d-1-a)-(d*(d-1))/(d-1-a)^3 - 2*d*(d-1)/(d-1-a)^4 + 3 * d * (d-1)^2 / (d-1-a)^5;
   % this case should not happen. But I am adding this to avoid NaN
end
end