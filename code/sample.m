function o = sample(plan, N)
% o = sample(plan, N)
%
%Input:
% plan - the sampling plan output by prepare.
% N - number of samples/steps
%
%Output:
% o - a structure continaing the following properties:
%   samples - a dim x N vector containing all the samples
%   samplesFullDim - a vector containing all the samples in the basis found
%   by the polytope class

%% Set default options
default.display = 0;
default.trajLength = 2;
default.minStepSize = 1e-4;
default.maxStepSize = 0.1;
default.maxRelativeStepSize = 0.2;
default.method = @implicitMidPoint;
default.rejectionSampling = false;
default.recordInterval = 1;
opts = setDefault(plan.opts,default); % add default if not specified

%% Sample
o = struct;
x = plan.initial;
o.samples = x;
ham = plan.ham; rejectedSamples = 0;
for i = 1:N
    if opts.display, fprintf('Iter %i\n', i); end
    rejected = 0;
    for trial = 1:10
        try
            z = ham.Generate(x);
            if opts.rejectionSampling, H1 = ham.H(z); end
            z = opts.method(z, ham, opts);
            if opts.rejectionSampling
                H2 = ham.H(z);
                if rand > min(1,exp(H1-H2))
                    z = NaN;
                    rejectedSamples = rejectedSamples + 1;
                    rejected = rejected + 1;
                end
            end
        catch
            z = NaN;
        end
        if ~isscalarnan(z), break; end
        if (rejected == 0)
            x = o.samples(:,max(1, size(o.samples,2)-trial));
            if opts.display, fprintf('Iter %i: ODE solver failed. Rewind %i steps\n', i, trial); end
        end
    end
    
    if isscalarnan(z)
        if opts.display
            if rejected == trial
                fprintf('Iter %i: rejection probability is too high. Set maxRelativeStepSize smaller.\n', i);
            end
        end
        fprintf('Solver ends prematurely.\n');
        break;
    end

    x = split2(z);
    
    % store the sample every few iterations
    if i > size(o.samples, 2) * opts.recordInterval
        o.samples = [o.samples zeros(size(o.samples))];
    end
    if mod(i,opts.recordInterval) == 0
        o.samples(:,i/opts.recordInterval) = x;
    end
end
o.samplesFullDim = o.samples(:,1:i/opts.recordInterval);
o.samples = plan.domain.T * o.samplesFullDim + plan.domain.y;