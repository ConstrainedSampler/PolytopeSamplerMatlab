function opts = default_options()
    opts = struct;
    opts.seed = 'shuffle';
    opts.nWorkers = 1;
    
    % Linear System Options
    opts.nSketch = 0;
    opts.solverThreshold = 1e-2;
    opts.extraHessian = 1e-20;
    opts.simdLen = 4;
    
    % ODE Options
    opts.odeMethod = @implicit_midpoint;
    opts.maxODEStep = 30;
    opts.implicitTol = 1e-5;
    
    % Termination Conditions
    opts.maxTime = 3600;
    opts.maxStep = +Inf;
	
    % HMC options
    opts.effectiveStepSize = 1;
    opts.initalStepSize = 0.2;
    opts.freezeMCMCAfterSamples = 100;
    opts.nRemoveInitialSamples = 5; % the number of samples we remove from the start
    
    % Module options
    opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'ProgressBar'};
    
    opts.DynamicStepSize = struct;
    opts.DynamicStepSize.maxConsecutiveBadStep = 10;
    opts.DynamicStepSize.targetODEStep = 10;
    opts.DynamicStepSize.shrinkFactor = 1.1;
    opts.DynamicStepSize.minStepSize = 0.001;
    opts.DynamicStepSize.warmUpStep = 10; % in terms of effective steps
    
    
    % We estimate the mixing time when the average accepted step 
    % = initialStep * stepMultiplier^k for k = 1, 2, ...
    opts.MixingTimeEstimator = struct;
    opts.MixingTimeEstimator.initialStep = 20; % in terms of effective steps
    opts.MixingTimeEstimator.stepMultiplier = 2;
    
    opts.MemoryStorage = struct;
    opts.MemoryStorage.recordsPerIndependentSample = 5;
    opts.MemoryStorage.minNumRecords = 100;
    
    % Presolve Options
    opts.presolve = struct;
    opts.presolve.runSimplify = true;
    opts.presolve.ipmMaxIter = 200;
    opts.presolve.ipmDualTol = 1e-12;
    opts.presolve.ipmDistanceTol = 1e-8; % we assume a constraint is tight if dist to constraint < distanceTol
    opts.presolve.splitDenseCols = 30;
    opts.presolve.removeFixedVariablesTol = 1e-12;
    
    % System Options
    opts.broadcastInterval = 0.5; % how often we sync between workers
    
    % Debug Options
    opts.rawOutput = false; % used only for debugging purpose, many functions for diagnostics does not work for raw output
    opts.profiling = false;
    opts.logging = []; % either a file name or a logging function of the form @(tag, msg, o) ...
end