function sample_test(debug, folders, problems, num_samples, func)
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
s.printFormat.mixing = '15f';
s.printFormat.pVal = '10f';
s.printFormat.preTime = '8.2f';
s.printFormat.stepSize = '10f';
s.printFormat.nStep = '10i';
s.printFormat.avgAcc = '15.3e';
s.testFunc = @(name) test_func(name, num_samples, debug, func);
s.test();
end

function o = test_func(name, num_samples, debug, func)

% load the problem and truncate it to make it bounded
P = loadProblem(name);
P_opts = default_options();
P_opts.maxTime = 3600*8;
P_opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DynamicWeight', 'DebugLogger'};
if debug
    P_opts.module{end+1} = 'ProgressBar';
end
d = size(P.lb,1);

switch func
   case 'uniform'
   case 'expoential'
      P.df = randn(d, 1);
   case 'normal'
      P.f = @(x) x'*x/2;
      P.df = @(x) x;
      P.ddf = @(x) ones(d,1);
end

sample_out = sample(P, num_samples, P_opts);

o = {};
o.m = size(sample_out.sampler.ham.A,1);
o.n = sample_out.sampler.ham.n;
o.nnz = nnz(sample_out.sampler.ham.A);
o.preTime = sample_out.prepareTime;
o.stepSize = sample_out.sampler.stepSize;
o.nStep = sample_out.totalStep;
o.avgAcc = mean(sample_out.averageAccuracy);
[o.pVal] = distribution_test(sample_out);
o.mixing = sample_out.sampler.mixingTime;

if (o.mixing < 500 && o.pVal > 0.005 && o.pVal < 0.995)
    o.success = 1;
else
    o.success = 0;
end
end

