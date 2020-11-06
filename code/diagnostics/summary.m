function smry = summary(spls)
%smry = summary(spls)
%compute summary of samples spls
%
%Input:
% spls - a dim x N vector, where N is the length of the chain.
%
%Output:
% smry - a table summarizing the mean and ess and rhat of spls.


st = std(spls, 0, 2);
m = mean(spls,2);
per25 =  prctile(spls, 25, 2);
per50 = prctile(spls, 50, 2);
per75 = prctile(spls, 75, 2);
ess = effective_sample_size(spls);
rh = rhat(spls);

[d, ~] = size(spls);

smry = table(m, st, per25, per50, per75, ess, rh, 'VariableNames',...
    {'mean','std','25%', '50%', '75%', 'n_ess', 'r_hat'}, ...
    'RowNames', 'samples['+string(1:d)+']');
end