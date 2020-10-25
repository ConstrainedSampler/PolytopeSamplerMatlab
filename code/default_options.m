function opts = default_options()
    opts = struct;
    opts.seed = 'shuffle';
    opts.maxTime = 3600;
    opts.maxStep = +Inf;
    opts.nSketch = 0;
    opts.maxODEStep = 20;
    opts.effectiveStepSize = 1;
    opts.initalStepSize = 0.2;
    opts.loggingFunc = @(tag, msg) 0;%fprintf('%s', msg);
    opts.freezeMCMCAfterSamples = 100;
    
    
    opts.module = {'MixingTimeEstimator', 'SampleStorage', 'DynamicRegularizer', 'DynamicStepSize', 'ProgressBar', 'DebugLogger'};
    
    opts.DynamicStepSize = struct;
    opts.DynamicStepSize.maxConsecutiveBadStep = 10;
    opts.DynamicStepSize.targetODEStep = 10;
    opts.DynamicStepSize.shrinkFactor = 1.1;
    opts.DynamicStepSize.minStepSize = 0.001;
    
    opts.MixingTimeEstimator = struct;
    opts.MixingTimeEstimator.startIter = 100;
    opts.MixingTimeEstimator.iterMultiplier = 2;
    
    opts.SampleStorage = struct;
    opts.SampleStorage.recordsPerIndependentSample = 5;
    opts.SampleStorage.minNumRecords = 100;
    opts.SampleStorage.rawOutput = false;
    
    
    
    % Other options
    opts.crudeSolverThreshold = 1e-6;
    opts.extraHessian = 1e-20;
    
    
    % ODE Options
    opts.odeMethod = @implicit_midpoint;
    opts.implicitTol = 1e-5;
    
    % Presolve Options
    opts.runSimplify = true;
    opts.ipmMaxIter = 100;
    opts.ipmDualTol = 1e-12;
    opts.ipmDistanceTol = 1e-8; % we assume a constraint is tight if dist to constraint < distanceTol
    opts.splitDenseCols = 30;
    opts.removeFixedVariablesTol = 1e-12;
end