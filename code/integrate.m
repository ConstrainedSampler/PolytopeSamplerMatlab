function W = integrate(problem, eps, opts, B, k)
%Input: a structure problem with the following fields
%  .Aineq
%  .bineq
%  .Aeq
%  .beq
%  .lb
%  .ub
%  .f
%  .df
%  .ddf
% describing a logconcave distribution given by
%   exp(-sum f_i(x_i))
%       over
%   {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}
% where f is given by a vector function of its first (df) and second (ddf)
% derivatives.
%
% Case 1: df is not defined
%   f(x) = 0.
% In this case, f, ddf must be empty. This is the uniform distribution. In
% this case the feasible region must be bounded.
%
% Case 2: df is a vector
%   f_i(x_i) = df_i x_i.
% In this case, f, ddf must be empty.
%
% Case 3: df is a function handle
%   f need to be defined as a function handle.
%   df need to be the derivative of f
%   ddf is optional. Providing ddf could improve the mixing time.
%
% eps - multiplicative error of the integral 
%
% opts - sampling options
%
% B - optional. Upper bound on the quantity max(f) - min(f)
%
% k - optional. Number of samples in each iteration. If k is provided, eps
% error is not guaranteed.
%
%Output:
% W - the integral of exp(-f) restricted on the polytope 
% {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}

%% Initialize parameters and compiling if needed
if (nargin <= 2)
    opts = default_options();
else
    opts = Setfield(default_options(), opts);
end

compile_solver(0); compile_solver(opts.simdLen);

%% Presolve
if isempty(opts.seed)
    opts.seed = randi(2^31);
end

if ischar(opts.logging) || isstring(opts.logging) % logging for Polytope
    fid = fopen(opts.logging, 'a');
    opts.presolve.logFunc = @(tag, msg) fprintf(fid, '%s', msg);
elseif ~isempty(opts.logging)
    opts.presolve.logFunc = opts.logging;
else
    opts.presolve.logFunc = @(tag, msg) 0;
end

polytope = Polytope(problem, opts);

if ischar(opts.logging) || isstring(opts.logging)
    fclose(fid);
end

d = polytope.n;
% fprintf('d is %d\n', d);
if (nargin <= 3)
    B = 2 * d;
end
num_of_iter = ceil(sqrt(d)*log(B));

if (nargin <= 4)
    k = 512/eps^2 * sqrt(d)*log(B);
end

W = 1;
a = 0;
f = polytope.f;
df = polytope.df;
ddf = polytope.ddf;

density = @(x) exp(-f(x));

problem.f = [];
problem.df = [];
problem.dff = [];

for i = 1:num_of_iter-1

    o = sample(problem, k, opts);
    W = W * mean(density(o.samples).^(1/B *(1+1/sqrt(d))^i - a));    
%     fprintf('W after iteration %d is %d\n', i, W);
    clear o
    a = 1/B * (1+1/sqrt(d))^i;
%     fprintf('a after iteration %d is %d\n', i, a);
    problem.f = @(x) a * f(x);
    problem.df = @(x) a * df(x);
    if isa(ddf, 'function_handle')
        problem.ddf = @(x) a * ddf(x);
    end
end


o = sample(problem, k, opts);
W = W * mean(density(o.samples).^(1 - a));
clear o

