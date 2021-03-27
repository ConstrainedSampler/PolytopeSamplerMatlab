function smry = summary(o)
%smry = summary(o)
%compute summary of samples
%
%Input:
% o - the samples object outputted by sample
%
%Output:
% smry - a table summarizing the mean and ess and rhat of samples.

% compute ess by summing over all chains
ess = 0;
for i = 1:numel(o.ess)
    ess_i = o.ess{i};
    for j = 1:size(ess_i,2)
        ess = ess + ess_i(:,j);
    end
end

% compute the rest of the summary
rh = rhat(o.chains);
st = std(o.samples, 0, 2);
m = mean(o.samples,2);
Y = prctile(o.samples, [25 50 75], 2);
per25 = Y(:, 1);
per50 = Y(:, 2);
per75 = Y(:, 3);
d = size(o.samples, 1);

smry = table(m, st, per25, per50, per75, ess, rh, 'VariableNames',...
    {'mean','std','first', 'second', 'third', 'ess', 'r_hat'}, ...
    'RowNames', 'samples['+string(1:d)+']');
end