function [samplePlan, flag] = prepare(problem, opts)
% [samplePlan, flag] = prepare(problem, opts)
% Output a sampling plan for the sample function to use
%
%Input:
% problem - a problem object describing the function to sample
% opts - a structure for options with the following properties
%  trajLength - trajectory length of each step in HMC (default:2)
%               Decrease this number to get better accuracy and robustness.
%  JLsize - number of direction used in estimating the drift (default:7)
%           Increase this number to get better accuracy.
%  maxRelativeStepSize - the maximum step size
%    (relative to the distance to the boundary) in the ODE solver (default:0.2)
%  method - the ODE solver (default:@implicitMidPoint)
%    Use @GaussLegendre4 for better accuracy
%  maxStepSize - the maximum step size in the ODE solver (default:0.1)
%  minStepSize - the minimum step size in the ODE solver  (default:1e-4)
%     If the ODE solver uses smaller step size multiple times,
%     the program quits
%  display - show the number of iterations (default:false)
%  recordInterval - record the sample every x iterations (default:1)
%  extraGaussian - restrict the step size by adding an x \cdot I term in the
%     metric. This can avoid the random walk runaway for unbounded.
%     Use it only for (nearly) unbounded body (default:0)
%Output:
% samplePlan - sampling plan for the sample function to use
% flag - 0 if successful
%      - 1 if the polytope is unbounded
%      - 2 if failed to compute the analytic center
%      - 3 if the polytope is unbounded and failed to compute the analytic center

%% Set default options
if ~exist('opts', 'var'), opts = struct; end
defaults.JLsize = 7;
defaults.extraGaussian = 0;
opts = setDefault(opts, defaults);

%% Prepare the polytope
P = Polytope(problem, opts);
P.simplify();
c = P.center;
flag = 0;

%% Detect Unbounded Polytope
w = estimateWidth(P);
if (max(w) > 1e9)
    warning('prepare:unbounded', 'Domain is unbounded. It may cause numerical problems. Make the body bounded or set extraGaussian = 1e-16 if necessary.');
    flag = flag + 1;
end

%% Create barrier using the df, ddf, dddf
phi = TwoSidedBarrier(P.lb, P.ub);
T = P.T; y = P.y;
T2 = T.^2; T3 = T.*T2;
if max(full(sum(T~=0,2))) > 1
    error('prepare:invalid_T', 'Each row of T in Polytope should contains at most 1 non-zero.');
end

if ~isempty(problem.f)
    if isa(problem.f,'function_handle')
        f = @(x) problem.f(T * x + y);
    else
        f = @(x) problem.f;
    end
    
    if any(size(f(c)) ~= [1 1])
        error('prepare:invalid_f', 'f should be a scalar function.');
    end
else
    f = @(x) 0;
end

if ~isempty(problem.df)
    if isa(problem.df,'function_handle')
        df = @(x) T'*problem.df(T * x + y);
    else
        tg = T' * problem.df;
        df = @(x) tg;
    end
    
    phi.extraGrad = df;
    if any(size(phi.extraGrad(c)) ~= [size(T,2) 1])
        error('prepare:invalid_df', 'df should be a n x 1 vector/function where n is the number of variables.');
    end
else
    zero_n = zeros(size(c));
    df = @(x) zero_n;
end

if ~isempty(problem.ddf)
    if isa(problem.ddf,'function_handle')
        phi.extraHess = @(x) T2'*problem.ddf(T * x + y);
    else
        th = T2' * problem.ddf;
        phi.extraHess = @(x) th;
    end
    
    if any(size(phi.extraHess(c)) ~= [size(T,2) 1])
        error('prepare:invalid_ddf', 'ddf should be a n x 1 vector/function where n is the number of variables.');
    end
end

if ~isempty(problem.dddf)
    if isa(problem.dddf,'function_handle')
        phi.extraTensor = @(x) T3'*problem.dddf(T * x + y);
    else
        th = T3' * problem.ddf;
        phi.extraTensor = @(x) th;
    end
    
    if any(size(phi.extraTensor(c)) ~= [size(T,2) 1])
        error('prepare:invalid_ddf', 'ddf should be a n x 1 vector/function where n is the number of variables.');
    end
end

%% Compute and check the analytic center
opts_ac = struct;
opts_ac.gaussianTerm = 1e-24;
phi.SetExtraHessian(opts.extraGaussian * ones(size(P.lb))); % this is important for free variables.
x = analyticCenter(P.A, P.b, phi, opts_ac);
Hinv = phi.HessianInv(x);
z = LinearSystemSolve(P.A, P.b-P.A*x, Hinv);
x = x + Hinv * (P.A' * z);
if ~(all(P.lb-1e-8<x) && all(x<1e-8+P.ub))
    warning('prepare:unbounded', 'Numerical problems occurred in finding an initial point. The body may be unbounded. Make the body bounded or set extraGaussian = 1e-16 if necessary.');
    flag = flag + 2;
end

samplePlan.domain = P;
samplePlan.ham = Hamiltonian(phi, P.A, P.b, f, df, opts.JLsize);
samplePlan.initial = x;
samplePlan.opts = opts;
end