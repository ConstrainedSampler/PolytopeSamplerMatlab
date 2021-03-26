function [ess] = effective_sample_size(x)
%ess = effective_sample_size(x)
%compute the effective sample sizes of each parameter (and each chain) in x
%
%Input:
% x - either a dim x N vector where N is the length of the chain,
%            a nChain x dim x N where nChain is the number of chains, or
%     a cell containing the matrices above
%
%Output:
% ess - either dim x 1 vector where ess(i) is the effective sample size of x(i,:)
%            a dim x nChain vector where ess(j,i) = ess(x(j,i,:)), or
%            a cell containing the outputs above

if iscell(x)
    % recurse over each element in the cell
    ess = cell(size(x));
    for i = 1:numel(x)
        ess{i} = effective_sample_size(x{i});
    end
else
    if numel(size(x)) == 3
        % recurse over each col in the tensor
        ess = zeros(size(x,2), size(x,1));
        for i = 1:size(x,1)
            ess(:,i) = effective_sample_size(squeeze(x(i,:,:)));
        end
    else
        N = size(x, 2); d = size(x, 1);
        Neven = N - mod(N,2);
        ess = zeros(d, 1);
        
        % compute ess over each row
        for i = 1:d
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
    end
end