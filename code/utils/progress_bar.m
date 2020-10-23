function progress_bar(total_samples, samples, time_spent, time_remain, acc_prob, step_size)
    persistent total_str_length
    if (total_samples == Inf)
        sample_num_length = 12;
        sample_num_field_length = 12;
    else
        sample_num_length = length(num2str(total_samples));
        sample_num_field_length = 2 * sample_num_length + 1;
    end
    bar_length = 25;
    
    if (nargin == 1)
        fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %', bar_length, 's | %', sample_num_field_length, 's | %8s | %8s\n');
        fprintf(fmt, 'Time spent', 'Time reamin', 'Progress', 'Samples', 'AccRate', 'StepSize');
        total_str_length = 0;
        return;
    end
    
    progress = samples / total_samples;
    
    if (total_str_length ~= 0)
        fprintf(repmat('\b',1, total_str_length));
    end
    
    if (total_samples == Inf)
        str = sprintf('%12s | %12s | %s | %12i | %8.6f | %8.6f', durationString(time_spent), durationString(time_remain), progressString(progress, bar_length), samples, acc_prob, step_size);
    else
        fmt = sprintf('%s%i%s%i%s', '%12s | %12s | %s | %', sample_num_length, 'i/%', sample_num_length, 'i | %8.6f | %8.6f');
        str = sprintf(fmt, durationString(time_spent), durationString(time_remain), progressString(progress, bar_length), samples, total_samples, acc_prob, step_size);
    end
    fprintf(str);
    total_str_length = length(str);
    
    function o = durationString(s)
        second = mod(floor(s),60);
        minute = mod(floor(s/60),60);
        hour = mod(floor(s/3600),24);
        day = floor(s/86400);
        o = sprintf('%02id:%02i:%02i:%02i', day, hour, minute, second);
    end

    function o = progressString(s, len)
        l = round(s*len);
        o = repmat('#', 1, l);
        o = [o, repmat(' ', 1, len - l)];
    end
end

