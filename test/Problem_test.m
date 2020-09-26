function Problem_test
s = TestSuite;
s.randomSeed = 1234;
s.nCores = 1;
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
s.randomSeed = 123456;
s.debug = 1;
s.test();
end

function o = test_func(name)
warning('off', 'Hamiltonian:Unbounded');

o = {};
P = loadProblem(name);

% Non-zero test vector to verify the OPT is same
if ~nonempty(P, 'df')
    P.df = randn(size(P.lb));
end

opts = optimoptions('linprog','Display','none');
P_opts = default_options();
P_opts.runSimplify = false;
P_opts.outputFunc = @(tag, msg, varargin) {};
P0 = Polytope(P, P_opts);
df = P0.df(zeros(P0.n,1));
x0 = linprog(P0.T' * df,[],[],P0.A,P0.b,P0.barrier.lb,P0.barrier.ub,opts);
if numel(x0) ~= 0
    o.opt = df' * (P0.T * x0 + P0.y);
else
    o.opt = NaN;
end

P_opts.runSimplify = true;
P1 = Polytope(P, P_opts);

o.feasible = P1.barrier.feasible(P1.center);
df = P1.df(zeros(P1.n,1));
x1 = linprog(P1.T' * df,[],[],P1.A,P1.b,P1.barrier.lb,P1.barrier.ub,opts);
if numel(x1) ~= 0
    o.optNew = df' * (P1.T * x1 + P1.y);
else
    o.optNew = NaN;
end

o.error = abs(o.opt-o.optNew)/(abs(o.opt)+1e-4);

% Problem 3
if contains(name,'netlib/')
    x = linprog(P.df, P.Aineq, P.bineq, P.Aeq, P.beq, P.lb, P.ub, struct('Display','none'));
    threshold = P.df' * x + abs(P.df)' * abs(x);
    P.Aineq = [P.Aineq; P.df'];
    P.bineq = [P.bineq; threshold];
    P.ub = min(P.ub, 2 * max(abs(x)));
    P.lb = max(P.lb, -2 * max(abs(x)));
end
P.df = [];
P3 = Polytope(P, P_opts);
c = P3.center;
o.minDist = min(c - P3.barrier.lb, P3.barrier.ub - c);
o.minDist = o.minDist ./ min(P3.barrier.ub-P3.barrier.lb, 1);
o.minDist = min(o.minDist);

o.m = size(P0.A,1); o.n = size(P0.A,2); o.nnz = nnz(P0.A);
o.mNew = size(P1.A,1); o.nNew = size(P1.A,2); o.nnzNew = nnz(P1.A);

if (o.error <= 1e-4 && o.feasible && o.minDist > 1e-5) % = is important for NaN case
    o.success = 1;
else
    o.success = 0;
end
end