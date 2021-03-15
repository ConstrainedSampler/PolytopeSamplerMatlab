initSampler

fid = fopen('demo.log', 'w');
P = loadProblem('basic/tv_ball@5');
P_opts = default_options();
P_opts.maxTime = 3600*8;
%P_opts.effectiveStepSize = 3;
%P_opts.initalStepSize = 0.1;
P_opts.SampleStorage.minNumRecords = 10000000;
P_opts.logFunc = @(tag, msg) fprintf(fid, '%s', msg); % Output the debug log to demo.log
%P_opts.module = {'MixingTimeEstimator', 'SampleStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DebugLogger', 'ProgressBar'};
P_opts.module = {'MixingTimeEstimator', 'SampleStorage', 'DebugLogger', 'ProgressBar'};
sample_out = sample(P, 100, P_opts);
fclose(fid);

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