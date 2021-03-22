function [rhat_val] = rhat(x)
%rhat_val = rhat(x)
%compute rhat of each parameter in the matrix x
%
%Input:
% x - a dim x N vector, where N is the length of the chain.
%
%Output:
% rhat - a dim x 1 vector where ess(i) is the effective sample size of x(i,:).

[~, N] = size(x);
if mod(N, 2) == 1
    x = x(:, 1:N-1);
    N = N-1;
end

N = N/2;
y = x(:, 1: N);
x = x(:, N+1: 2*N);

b_div_n = var([mean(x, 2), mean(y, 2)], 0, 2);
w = mean([var(x, 0, 2), var(y, 0, 2)], 2) + eps;

sig_2p = (N-1)/N .* w + b_div_n;

rhat_val = 3/2 * sig_2p ./w - (N-1)/ (2*N);
end