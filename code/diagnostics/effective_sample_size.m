function [ess] = effective_sample_size(x)
%ess = effective_sample_size(x)
%compute the effective sample sizes of each parameter in the matrix x
%
%Input:
% x - a dim x N vector, where N is the length of the chain.
%
%Output:
% ess - a dim x 1 vector where ess(i) is the effective sample size of x(i,:).

n = size(x, 2);
n_even = n - mod(n,2);
ess = zeros(size(x,1), 1);

for i = 1:size(x,1)
    % normalize i-th row
    x_ = x(i,:);
    
    m = mean(x_); 
    x_ = x_ - m;
    
    v = mean(x_.^2) + 1e-16 * abs(m).^2; % Avoid dividing by 0
    x_ = x_ / sqrt(v);
    
    % compute autocorrelation via Wiener–Khinchin theorem
    r = ifft(abs(fft(x_,2*n)).^2); % power spectral density
    ac = real(r(1:n_even))/n; % ess formula assume vector length is even
    
    % Geyer's monotone estimator
    minAC = cummin(ac(:, 1:2:end) + ac(:, 2:2:end));
    ess(i) = n/max(1,2*sum(minAC.*(minAC>0)) -1);
end
end