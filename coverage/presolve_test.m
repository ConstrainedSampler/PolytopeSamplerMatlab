function presolve_test(debug, folders, problems)
s = TestSuite;
if nargin >= 2 && ~isempty(folders)
    s.problemFilter.folders = folders;
end
if nargin >= 3 && ~isempty(problems)
    s.problems = problems;
end
s.randomSeed = 123456;
s.nCores = +Inf;
s.debug = debug;
s.printFormat.m = '8i';
s.printFormat.n = '8i';
s.printFormat.nnz = '10i';
s.printFormat.mNew = '8i';
s.printFormat.nNew = '8i';
s.printFormat.nnzNew = '10i';
s.printFormat.opt = '14e';
s.printFormat.optNew = '14e';
s.printFormat.error = '10.3e';
s.printFormat.feasible = '10i';
s.printFormat.minDist = '10.3e';
s.testFunc = @test_func;
s.test();
end

function o = test_func(name)
o = {};
P = loadProblem(name);
rng(123456); % make sure it is reproducible

%% Test 1: Check if the solution remains the same.

if isempty(P.df)
    P.df = randn(size(P.lb));
end

% first solve it using matlab LP solver
P_opts = default_options();
P_opts.presolve.runSimplify = false;
P.center = P.df; % avoid Polytope raise error of not finding a center
P_opts.presolve.logFunc = @(tag, msg) 0;
P0 = Polytope(P, P_opts);
P.center = [];
df = P0.df;

opts = optimoptions('linprog','Display','none');
x0 = linprog(P0.T' * df, [], [], P0.A, P0.b, P0.barrier.lb, P0.barrier.ub, opts);
if numel(x0) ~= 0
    o.opt = df' * (P0.T * x0 + P0.y);
else
    o.opt = NaN;
end

% simplify the polytope and solve it again
P_opts.presolve.runSimplify = true;
P1 = Polytope(P, P_opts);
df = P1.df;

x1 = linprog(P1.T' * df, [], [], P1.A, P1.b, P1.barrier.lb, P1.barrier.ub, opts);
if numel(x1) ~= 0
    o.optNew = df' * (P1.T * x1 + P1.y);
else
    o.optNew = NaN;
end

o.feasible = P1.barrier.feasible(P1.center);
o.error = abs(o.opt-o.optNew)/(abs(o.opt)+1e-4);

%% Test 2: check if the analytic center is deep inside
P.df = zeros(size(P.lb));
P3 = Polytope(P, P_opts);
c = P3.center;
o.minDist = min(c - P3.barrier.lb, P3.barrier.ub - c);
o.minDist = o.minDist ./ min(P3.barrier.ub-P3.barrier.lb, 1);
o.minDist = min(o.minDist);

o.m = size(P0.A,1); o.n = size(P0.A,2); o.nnz = nnz(P0.A);
o.mNew = size(P1.A,1); o.nNew = size(P1.A,2); o.nnzNew = nnz(P1.A);

if (o.error <= 1e-4 && o.feasible && o.minDist > 1e-6) % = is important for NaN case
    o.success = 1;
else
    o.success = 0;
end
end