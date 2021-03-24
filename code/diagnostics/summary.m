function smry = summary(o)
%smry = summary(o)
%compute summary of samples
%
%Input:
% o - output of sample
%
%Output:
% smry - a table summarizing the mean and ess and rhat of spls.

% compute ess by summing over all chains
ess = 0;
for i = 1:numel(o.ess)
    ess_i = o.ess{i};
    for j = 1:size(ess_i,2)
        ess = ess + ess_i(:,j);
    end
end

% compute the rest of the info
spls = o.samples;
rh = rhat(spls);

merged_spls = cell2mat(spls);
st = std(merged_spls, 0, 2);
m = mean(merged_spls,2);
Y = prctile(merged_spls, [25 50 75], 2);
per25 = Y(:, 1);
per50 = Y(:, 2);
per75 = Y(:, 3);
d = size(merged_spls, 1);
clear merged_spls

smry = table(m, st, per25, per50, per75, ess, rh, 'VariableNames',...
    {'mean','std','25%', '50%', '75%', 'ess', 'r_hat'}, ...
    'RowNames', 'samples['+string(1:d)+']');
end