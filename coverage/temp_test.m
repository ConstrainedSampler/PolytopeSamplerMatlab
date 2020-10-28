initSampler


P = loadProblem('basic/box@10000');
P_opts = default_options();
P_opts.maxTime = 3600*8;
P_opts.SampleStorage.minNumRecords = 100000000;
P_opts.module = {'MixingTimeEstimator', 'SampleStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DebugLogger'};
sample_out = sample(P, 100, P_opts);

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