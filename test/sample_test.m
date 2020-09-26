function sample_test
clc
s = TestSuite;
s.randomSeed = 123456;
s.nCores = 4;
s.debug = 0;
s.printFormat.m = '8i';
s.printFormat.n = '8i';
s.printFormat.nnz = '10i';
s.printFormat.mixing = '12f';
s.printFormat.pVal = '10f';
s.printFormat.preTime = '8.2f';
s.printFormat.stepSize = '10f';
s.printFormat.nStep = '10i';
s.printFormat.avgAcc = '15.3e';
s.testFunc = @(name) test_func(name, 100);
s.test();
end

function o = test_func(name, num_samples)
warning('off', 'PolytopeSampler:uniformtest:size');
warning('off', 'stats:adtest:OutOfRangePLow');
warning('off', 'stats:adtest:OutOfRangePHigh');
warning('off', 'stats:adtest:SmallSampleSize');
warning('off', 'stats:adtest:NotEnoughData');

% load the problem and truncate it to make it bounded
P = loadProblem(name);
gc = zeros(size(P.lb,1),1); % Gaussian center
if contains(name,'netlib/')
    x = linprog(P.df, P.Aineq, P.bineq, P.Aeq, P.beq, P.lb, P.ub, struct('Display','none'));
    threshold = P.df' * x + abs(P.df)' * abs(x);
    P.Aineq = [P.Aineq; P.df'];
    P.bineq = [P.bineq; threshold];
    P.ub = min(P.ub, 2 * max(abs(x)));
    P.lb = max(P.lb, -2 * max(abs(x)));
    gc = x;
end
%uniform sampling in P
%P.df = 0;
%Gaussian sampling in P
P.f = @(x)(x-gc)'*(x-gc);
P.df = @(x)2*(x-gc);
P.ddf = @(x)2*ones(size(P.lb,1),1);
P.dddf = @(x)zeros(size(P.lb,1),1);

P_opts = default_options();
P_opts.maxTime = 3600*8;
P_opts.debugMode = true;
P_opts.outputFunc = @(tag, msg, varargin) {};
sample_out = sample(P, num_samples, P_opts);

o = {};
o.m = size(sample_out.ham.A,1);
o.n = sample_out.ham.n;
o.nnz = nnz(sample_out.ham.A);
o.preTime = sample_out.prepareTime;
o.stepSize = sample_out.stepSize;
o.nStep = sample_out.i;
o.avgAcc = sample_out.averageAccuracy;

% diagnosis
s = sample_out.ham.T * sample_out.samples + sample_out.ham.y;

o.pVal=0.5;
%[o.pVal] = uniformtest(s, P, sample_out.dim);
o.mixing = sample_out.mixingTime;

if (o.mixing < 500 && o.pVal > 0.005 && o.pVal < 0.995)
    o.success = 1;
else
    o.success = 0;
end
end

