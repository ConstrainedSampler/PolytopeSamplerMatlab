function smry = summary(o)
%smry = summary(o)
%compute summary of samples
%
%Input:
% o - the samples object outputted by sample
%
%Output:
% smry - a table summarizing the mean and rhat of samples.

% compute ess by summing over all chains

% compute the rest of the summary
rh = rhat(o.samples);
st = std(o.samples, 0, 2);
m = mean(o.samples,2);
d = size(o.samples, 1);

if isempty(ver('stats'))
    smry = table(m, st, rh, 'VariableNames',...
        {'mean', 'std', 'r_hat'}, ...
        'RowNames', 'samples['+string(1:d)+']');
else
    Y = prctile(o.samples, [25 50 75], 2);
    per25 = Y(:, 1);
    per50 = Y(:, 2);
    per75 = Y(:, 3);

    smry = table(m, st, per25, per50, per75, rh, 'VariableNames',...
        {'mean', 'std', 'Q1', 'Q2', 'Q3', 'r_hat'}, ...
        'RowNames', 'samples['+string(1:d)+']');
end
end