function uniform_sample_test(debug, folder, problems)
s = TestSuite;
if nargin >= 2 && ~isempty(folder)
    s.problemFilter.folder = folder;
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
s.printFormat.mixing = '12f';
s.printFormat.pVal = '10f';
s.printFormat.preTime = '8.2f';
s.printFormat.stepSize = '10f';
s.printFormat.nStep = '10i';
s.printFormat.avgAcc = '15.3e';
s.testFunc = @(name) test_func(name, 200);
s.test();
end

function o = test_func(name, num_samples)

% load the problem and truncate it to make it bounded
P = loadProblem(name);
P_opts = default_options();
P_opts.maxTime = 3600*8;
P_opts.module = {'MixingTimeEstimator', 'SampleStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DebugLogger'};
sample_out = sample(P, num_samples, P_opts);

o = {};
o.m = size(sample_out.sampler.ham.A,1);
o.n = sample_out.sampler.ham.n;
o.nnz = nnz(sample_out.sampler.ham.A);
o.preTime = sample_out.prepareTime;
o.stepSize = sample_out.sampler.stepSize;
o.nStep = sample_out.totalStep;
o.avgAcc = sample_out.averageLinearSystemAccuracy;
[o.pVal] = uniformtest(sample_out);
o.mixing = sample_out.sampler.mixingTime;

if (o.mixing < 500 && o.pVal > 0.005 && o.pVal < 0.995)
    o.success = 1;
else
    o.success = 0;
end
end

