function opts = default_options()
%Defaults set to no JL, and no weightedBarrier
    opts = struct;
    opts.seed = 'shuffle';
    opts.nSketch = 0;
    opts.adaptiveStepSize = true;
    opts.weightedBarrier = false;
    opts.dynamicBound = true;
    opts.crudeSolverThreshold = 1e-6;
    opts.extraHessian = 1e-20;
    
    
    % Algorithm Options
    opts.effectiveStepSize = 1;
    opts.shrinkFactor = 1.1;
    opts.initalStepSize = 0.2;
    
    % ODE Options
    opts.odeMethod = @implicit_midpoint;
    opts.maxODEStep = 20;
    opts.targetODEStep = 10;
    opts.implicitTol = 1e-5;
    
    % Presolve Options
    opts.runSimplify = true;
    opts.ipmMaxIter = 100;
    opts.ipmDualTol = 1e-12;
    opts.ipmDistanceTol = 1e-8; % we assume a constraint is tight if dist to constraint < distanceTol
    opts.splitDenseCols = 30;
    opts.removeFixedVariablesTol = 1e-12;
    
    % Output Options
    opts.recordsPerIndependentSample = 10;
    opts.outputFunc = @(tag, msg, varargin) fprintf(msg, varargin{:});
    opts.rawOutput = false;
    
    % Terminiation Condition
    opts.checkESSIter = 100;
    opts.maxTime = 3600;
    opts.maxStep = +Inf;
    opts.minStepSize = 0.001;
    opts.minAcceptedSteps = 100;
    
    % Restart/Shrink Condition
    opts.consecutiveBadStepLimit = 10;
    opts.freezeMCMCAfterSamples = 100;
end