classdef ProgressBar < handle
    properties
        sampler
        opts
        
        startTime
        acceptedStep = 0;
        nextBackspaceLength
        barLength = 25
        nSamplesTextLength = 12
        nSamplesFieldLength = 12
        nBackspacePerPrint
        
        refreshInterval = 0.2
        lastRefresh
    end
    
    methods
        function o = ProgressBar(sampler)
            o.sampler = sampler;
            o.startTime = tic;
            o.lastRefresh = tic;
            o.acceptedStep = 0;
            if (sampler.N ~= Inf)
                o.nSamplesTextLength = max(length(num2str(sampler.N)*2),3);
                o.nSamplesFieldLength = 2 * o.nSamplesTextLength + 1;
            end
            
            if sampler.labindex == 1
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %', o.barLength, 's | %', o.nSamplesFieldLength, 's | %8s | %8s | %8s\n');
                s = sprintf(fmt, 'Time spent', 'Time reamin', 'Progress', 'Samples', 'AccRate', 'StepSize', 'MixTime');
                disp(s);
                
                if (sampler.nWorkers > 1)
                    o.nBackspacePerPrint = 3;
                else
                    o.nBackspacePerPrint = 1;
                end
                o.nextBackspaceLength = o.nBackspacePerPrint;
            end
        end
        
        function o = initialize(o)
            
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            o.acceptedStep = o.acceptedStep + sum(o.sampler.accept);
            if toc(o.lastRefresh) < o.refreshInterval, return, end
            if o.sampler.labindex == 1
                o.refresh_bar();
            end
        end
        
        function o = refresh_bar(o)
            s = o.sampler; k = size(s.samples,1);
            o.lastRefresh = tic;
            prob = (o.acceptedStep/k) / s.i;
            timeSpent = toc(o.startTime);
            avgMixingTime = (size(s.samples,1) * s.nWorkers) / s.sampleRate;
            if isnan(s.sampleRate)
                nSamples = 0;
                timeRemain = Inf;
            else
                nSamples = floor(max(s.i * s.sampleRate, s.totalNumSamples));
                progress = min(1, nSamples / s.N);
                timeRemain = timeSpent / progress * (1-progress);
            end
            timeRemain = min(timeRemain, s.opts.maxTime - timeSpent);
            timeRemain = min(timeRemain, (s.opts.maxStep - s.i) / s.i * timeSpent);
            progress = timeSpent/(timeRemain+timeSpent);
            
            str1 = sprintf(repmat('\b', 1, o.nextBackspaceLength));
            
            if (o.sampler.N == Inf)
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesFieldLength, 'i | %8.6f | %8.6f | %8.1f');
                str2 = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, prob, s.stepSize, avgMixingTime);
            else
                fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', o.nSamplesTextLength, 'i/%', o.nSamplesTextLength, 'i | %8.6f | %8.6f | %8.1f');
                str2 = sprintf(fmt, durationString(timeSpent), durationString(timeRemain), progressString(progress, o.barLength), nSamples, s.N, prob, s.stepSize, avgMixingTime);
            end
            
            disp([str1 str2]);
            o.nextBackspaceLength = length(str2)+o.nBackspacePerPrint;

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
        
        function o = sync(o)
            
        end
        
        function o = finalize(o)
            if o.sampler.labindex == 1
                o.refresh_bar();
                fprintf('\nDone!\n');
            end
        end
    end
end