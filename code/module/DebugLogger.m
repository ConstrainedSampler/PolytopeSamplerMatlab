classdef DebugLogger < handle
    % Module for output the following information in the sample output
    % - acceptedStep, the number of accepted step for each chain
    % - totalStep, the total number of step taken
    % - numCholesky(1:end-1), the number of high precision cholesky decompositions we performed
    % - numCholesky(end), the number of cholesky decompositions (in any precision) we performed
    % - sampler, the sampler object we used
    %
    % If the sampler fails, dump all variables to 'dump_#_ignore.mat'
    methods
        function o = DebugLogger(s)
            s.output.averageAccuracy = 0;
        end
        
        function o = step(o, s)
            s.output.averageAccuracy = s.output.averageAccuracy + s.ham.solver.accuracy;
            if (any(s.ham.solver.accuracy > s.opts.solverThreshold))
                s.log('Solver:inaccurate', 'low double accuracy %s.\n', num2str(s.ham.solver.accuracy'));
            end
        end
        
        function o = finalize(o, s)
            s.output.acceptedStep = s.acceptedStep;
            s.output.totalStep = s.i; % ignore the warm up phase
            s.output.numCholesky = s.ham.solver.getDecomposeCount();
            s.output.averageAccuracy = s.output.averageAccuracy / s.i;
            s.output.sampler = s;
            if s.terminate == 3
                save(sprintf('dump_%i_ignore.mat', s.nWorkers));
            end
        end
    end
end