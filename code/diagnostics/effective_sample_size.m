function [ess] = effective_sample_size(x)
%ess = effective_sample_size(x)
%compute the effective sample sizes of each parameter in the matrix x
%
%Input:
% x - a dim x N vector, where N is the length of the chain.
%
%Output:
% ess - a dim x 1 vector where ess(i) is the effective sample size of x(i,:).

s = size(x);
N = s(end);
Neven = N - mod(N,2);
l = prod(s(1:end-1));
ess = zeros(l, 1);
x = reshape(x, [l, N]);

for i = 1:l
    % normalize i-th row
    x_ = x(i,:);
    m = mean(x_); 
    x_ = x_ - m;
    
    v = mean(x_.^2) + 1e-16 * abs(m).^2; % Avoid dividing by 0
    x_ = x_ / sqrt(v);
    
    % compute autocorrelation via Wiener–Khinchin theorem
    r = ifft(abs(fft(x_,2*N)).^2); % power spectral density
    ac = real(r(1:Neven))/N; % ess formula assume vector length is even
    
    % Geyer's monotone estimator
    minAC = cummin(ac(:, 1:2:end) + ac(:, 2:2:end));
    ess(i) = N/max(1,2*sum(minAC.*(minAC>0)) -1);
end

ess = reshape(ess, [s(1:end-1) 1]);
end