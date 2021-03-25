initSampler_ignore

%profile on;
P = loadProblem('metabolic/Recon1');
%P = loadProblem('netlib/agg');
opts = default_options();
opts.seed = 1;
opts.simdLen = 4;
opts.outputFormat = 'separate';
opts.nWorkers = 1;
opts.maxTime = 300;
%opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicStepSize', 'ProgressBar'};
%opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'ProgressBar'};
%opts.maxStep = 4000;
% requirement, it will simply need to non-convergence

opts.module{end+1} = 'DebugLogger';
%opts.initalStepSize = 0.05;
opts.logging = 'ttt_ignore.log';
o = sample(P, 1000, opts);
%[pVal] = uniformtest(o, struct('toPlot', true));
%profile report;

% old version
%   Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
% 00d:00:02:44 | 00d:00:00:00 | ######################### | 105/100 | 0.828560 | 0.200000 |    360.1
% 0.0043386 per step

% new version
%  Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
%00d:00:02:25 | 00d:00:00:00 | ######################### | 105/100 | 0.828037 | 0.200000 |    308.7
% 0.0044734 per step 
% 0.0042340 per step (after limiting use of T) (For recon1, only 50% time spend on PackedChol) 

% matlab version
%  Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
%00d:00:04:54 | 00d:00:00:00 | ######################### | 110/100 | 0.815321 | 0.200000 |    377.9
% 0.0070725 per step

% fake simd version
%  Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
%00d:00:12:22 | 00d:00:00:00 | ######################### | 421/400 | 0.906054 | 0.165289 |    278.4
% 0.0063307 per step

% real simd version
%  Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
%00d:00:03:10 | 00d:00:00:00 | ######################### | 113/100 | 0.886729 | 0.181818 |    385.5
% 0.0043616 per step
% 0.0040166 per step (without profile)

%% Recon3 (limit 5 min)
% 00d:00:04:24<----time bug % i = 1904
% 00d:00:03:35<----time bug % i = 1904