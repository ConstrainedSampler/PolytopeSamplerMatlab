classdef ProgressBar < handle
    % Module for printing out progress bar on Matlab console
    
    properties
        barLength = 25 % the length of the field "Progress"
        nSamplesTextLength = 12 % the length of "samples number"
        nSamplesFieldLength = 12 % the length of the field "Est Samples"
        nextBackspaceLength
        nExtraBackspacePerPrint
        
        refreshInterval = 0.2
        lastRefresh = tic();
    end
    
    methods
        function o = ProgressBar(s)
            if (s.N ~= Inf)
                o.nSamplesTextLength = max(length(num2str(s.N)), 5);
                o.nSamplesFieldLength = 2 * o.nSamplesTextLength + 1;
            end
            
            if s.labindex == 1
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %', o.barLength, 's | %', o.nSamplesFieldLength, 's | %8s | %8s | %8s\n');
                s = sprintf(fmt, 'Time spent', 'Time reamin', 'Progress', 'Est Samples', 'AccRate', 'StepSize', 'MixTime');
                disp(s);
                
                if (s.nWorkers > 1)
                    o.nExtraBackspacePerPrint = 3;
                else
                    o.nExtraBackspacePerPrint = 1;
                end
                o.nextBackspaceLength = o.nExtraBackspacePerPrint;
            end
        end
        
        function o = step(o, s)
            if toc(o.lastRefresh) < o.refreshInterval, return, end
            if s.labindex == 1
                o.refresh_bar();
            end
        end
        
        function o = finalize(o, s)
            if s.labindex == 1
                o.refresh_bar();
                fprintf('Done!\n');
            end
        end
        
        function o = refresh_bar(o, s)
            prob = mean(s.acceptedStep) / s.i;
            timeSpent = toc(s.startTime);                                    % s.startTime = time after presolve
            
            if isnan(s.sampleRate)
                avgMixingTime = NaN;
                nSamples = 0;
                timeRemain = Inf;
            else
                avgMixingTime = (size(s.x,1) * s.nWorkers) / s.sampleRate;
                nSamples = floor(max(s.i * s.sampleRate, s.totalNumSamples));
                r = min(1, nSamples / s.N);
                timeRemain = (timeSpent / r) * (1-r);
            end
            timeRemain = min(timeRemain, s.opts.maxTime - s.opts.startTime); % s.opts.startTime = time before presolve
            timeRemain = min(timeRemain, (s.opts.maxStep - s.i) / s.i * timeSpent);
            progress = timeSpent/(timeRemain+timeSpent);
            
            str1 = sprintf(repmat('\b', 1, o.nextBackspaceLength));
            
            if (s.N == Inf)
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesFieldLength, 'i | %8.6f | %8.6f | %8.1f');
                str2 = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, prob, s.stepSize, avgMixingTime);
            else
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesTextLength, 'i/%', o.nSamplesTextLength, 'i | %8.6f | %8.6f | %8.1f');
                str2 = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, s.N, prob, s.stepSize, avgMixingTime);
            end
            
            disp([str1 str2]);
            o.nextBackspaceLength = length(str2) + o.nExtraBackspacePerPrint;
            o.lastRefresh = tic;
            
            function o = durationString(str)
                if isnan(str) || isinf(str)
                    o = '         NaN';
                else
                    str = max(str,0);
                    second = mod(floor(str),60);
                    minute = mod(floor(str/60),60);
                    hour = mod(floor(str/3600),24);
                    day = floor(str/86400);
                    o = sprintf('%02id:%02i:%02i:%02i', day, hour, minute, second);
                end
            end

            function o = progressString(str, len)
                str = max(min(str, 1),0);
                l = round(str*len);
                o = repmat('#', 1, l);
                o = [o, repmat(' ', 1, len - l)];
            end
        end
    end
end