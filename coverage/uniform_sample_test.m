function uniform_sample_test
s = TestSuite;
s.randomSeed = 123456;
s.nCores = +Inf;
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
s.testFunc = @(name) test_func(name, 200);
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
P_opts = default_options();
P_opts.maxTime = 3600*8;
P_opts.outputFunc = @(tag, msg, varargin) {};
sample_out = sample(P, num_samples, P_opts);

o = {};
o.m = size(sample_out.ham.A,1);
o.n = sample_out.ham.n;
o.nnz = nnz(sample_out.ham.A);
o.preTime = sample_out.prepareTime;
o.stepSize = sample_out.stepSize;
o.nStep = sample_out.i;
o.avgAcc = sample_out.averageLinearSystemAccuracy;
[o.pVal] = uniformtest(sample_out);
o.mixing = sample_out.mixingIter;

if (o.mixing < 500 && o.pVal > 0.005 && o.pVal < 0.995)
    o.success = 1;
else
    o.success = 0;
end
end

