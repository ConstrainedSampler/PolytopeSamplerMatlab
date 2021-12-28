classdef DynamicStepSize < handle
    % Module for dynamically choosing the step size
    properties
        opts
        
        consecutiveBadStep
        iterSinceShrink
        rejectSinceShrink
        ODEStepSinceShrink
        effectiveStep
        warmupFinished = false
    end
    
    methods
        function o = DynamicStepSize(s)
            o.opts = s.opts.DynamicStepSize;
            o.consecutiveBadStep = 0;
            o.iterSinceShrink = 0;
            o.rejectSinceShrink = 0;
            o.ODEStepSinceShrink = 0;
            o.effectiveStep = 0;
            
            if o.opts.warmUpStep > 0
                s.stepSize = 1e-3;
            else
                o.warmupFinished = true;
            end
        end
        
        function o = step(o, s)
            % Warmup phase
            bad_step = s.prob < 0.5 | s.ODEStep == s.opts.maxODEStep;
            o.consecutiveBadStep = bad_step .* o.consecutiveBadStep + bad_step;
            warmupRatio = mean(s.nEffectiveStep) / o.opts.warmUpStep;
            
            if warmupRatio < 1 && ~o.warmupFinished && max(o.consecutiveBadStep) < o.opts.maxConsecutiveBadStep
                s.stepSize = s.opts.initalStepSize * min(warmupRatio+1e-2, 1);
                s.momentum = 1 - min(1, s.stepSize / s.opts.effectiveStepSize);
                return;
            end
            
            if (~o.warmupFinished)
                s.i = 1;
                s.acceptedStep = 0;
                s.nEffectiveStep = 0;
                s.chains = zeros(s.opts.simdLen, s.ham.n, 0);
                o.warmupFinished = true;
            end
            
            o.iterSinceShrink = o.iterSinceShrink + 1;
            o.rejectSinceShrink = o.rejectSinceShrink + 1-s.prob;
            o.ODEStepSinceShrink = o.ODEStepSinceShrink + s.ODEStep;
            
            % Shrink the step during the warmup phase
            if (~s.freezed)
                shrink = 0;
                shiftedIter = o.iterSinceShrink + 20 / (1-s.momentum);
                
                targetProbability = (1-s.momentum)^(2/3)/4;
                if (max(o.rejectSinceShrink) > targetProbability  * shiftedIter)
                    shrink = sprintf('Failure Probability is %.4f, which is larger than the target %.4f', max(o.rejectSinceShrink) / o.iterSinceShrink, targetProbability);
                end
                
                if (max(o.consecutiveBadStep) > o.opts.maxConsecutiveBadStep)
                    shrink = sprintf('Consecutive %i Bad Steps', max(o.consecutiveBadStep));
                end
                
                if (max(o.ODEStepSinceShrink) > o.opts.targetODEStep * shiftedIter)
                    shrink = sprintf('ODE solver requires %.4f steps in average, which is larger than the target %.4f', max(o.ODEStepSinceShrink) / o.iterSinceShrink, o.opts.targetODEStep);
                end
                
                if ischar(shrink)
                    o.iterSinceShrink = 0;
                    o.rejectSinceShrink = 0;
                    o.ODEStepSinceShrink = 0;
                    o.consecutiveBadStep = 0;
                    
                    s.stepSize = s.stepSize / o.opts.shrinkFactor;
                    s.momentum = 1 - min(0.999, s.stepSize / s.opts.effectiveStepSize);
                    
                    s.log('DynamicStepSize:step', 'Step shrinks to h = %f due to %s\n', s.stepSize, shrink);
                    
                    if s.stepSize < o.opts.minStepSize
                        s.log('warning', 'Algorithm fails to converge even with step size h = %f.\n', s.stepSize);
                        s.terminate = 3;
                    end
                end
                
                o.iterSinceShrink = o.iterSinceShrink + 1;
            elseif max(o.consecutiveBadStep) > o.opts.maxConsecutiveBadStep
                s.x = ones(size(s.x,1),1) * mean(s.chains, [1 3]);
                s.v = s.ham.resample(s.x, zeros(size(s.x)));
                s.log('DynamicStepSize:step', 'Sampler reset to the center gravity due to consecutive bad steps.\n');
                
                o.iterSinceShrink = 0;
                o.rejectSinceShrink = 0;
                o.ODEStepSinceShrink = 0;
                o.consecutiveBadStep = 0;
            end
        end
    end
end