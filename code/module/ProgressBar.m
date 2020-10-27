classdef ProgressBar < handle
    properties
        sampler
        opts
        
        startTime
        acceptedStep = 0;
        lastDisplayLength = 0
        barLength = 25
        nSamplesTextLength = 12
        nSamplesFieldLength = 12
        
        refreshInterval = 0.2
        lastRefresh
    end
    
    methods
        function o = ProgressBar(sampler)
            o.sampler = sampler;
            o.startTime = tic;
            o.lastRefresh = tic;
            sampler.acceptedStep = 0;
            if (sampler.N ~= Inf)
                o.nSamplesTextLength = max(length(num2str(sampler.N)*2),3);
                o.nSamplesFieldLength = 2 * o.nSamplesTextLength + 1;
            end
            
            fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %', o.barLength, 's | %', o.nSamplesFieldLength, 's | %8s | %8s\n');
            fprintf(fmt, 'Time spent', 'Time reamin', 'Progress', 'Samples', 'AccRate', 'StepSize');
        end
        
        function o = initialize(o)
            
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            o.acceptedStep = o.acceptedStep + o.sampler.accept;
            if toc(o.lastRefresh) < o.refreshInterval, return, end
            o.refresh_bar();
        end
        
        function o = refresh_bar(o)
            s = o.sampler;
            o.lastRefresh = tic;
            prob = o.acceptedStep / s.i;
            timeSpent = toc(o.startTime);
            if isnan(s.mixingTime) || s.mixingTime <= 0.5
                progress = 0;
                nSamples = 0;
                timeRemain = Inf;
            else
                nSamples = floor(s.i / s.mixingTime);
                progress = min(1, nSamples / s.N);
                timeRemain = timeSpent / progress * (1-progress);
            end
            timeRemain = min(timeRemain, s.opts.maxTime - timeSpent);
            timeRemain = min(timeRemain, (s.opts.maxStep - s.i) / s.i * timeSpent);
            progress = timeSpent/(timeRemain+timeSpent);
            
            if (o.lastDisplayLength ~= 0)
                fprintf(repmat('\b', 1, o.lastDisplayLength));
            end
            
            if (o.sampler.N == Inf)
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesFieldLength, 'i | %8.6f | %8.6f');
                str = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, prob, s.stepSize);
            else
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesTextLength, 'i/%', o.nSamplesTextLength, 'i | %8.6f | %8.6f');
                str = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, s.N, prob, s.stepSize);
            end
            
            fprintf(str);
            o.lastDisplayLength = length(str);

            function o = durationString(s)
                if isnan(s) || isinf(s)
                    o = '         NaN';
                else
                    s = max(s,0);
                    second = mod(floor(s),60);
                    minute = mod(floor(s/60),60);
                    hour = mod(floor(s/3600),24);
                    day = floor(s/86400);
                    o = sprintf('%02id:%02i:%02i:%02i', day, hour, minute, second);
                end
            end

            function o = progressString(s, len)
                s = max(min(s, 1),0);
                l = round(s*len);
                o = repmat('#', 1, l);
                o = [o, repmat(' ', 1, len - l)];
            end
        end
        
        function o = finalize(o)
            o.refresh_bar();
            fprintf('\nDone!\n');
        end
    end
end