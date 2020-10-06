initSampler

%% Example 1: Sample uniform from a simplex
P = struct; d = 100;
P.Aineq = ones(1, d);
P.bineq = 1;
P.lb = zeros(d, 1);

o = sample(P, 100); % Number of samples = 100
s = thin_samples(o.samples);  % extract "independent" samples from the dependent samples
histogram(sum(s), 0.9:0.005:1)
title('distribution of l1 norm of simplex');
[pVal] = uniformtest(o, struct('toPlot', true));
drawnow()

%% Example 2: Sample uniform from Birkhoff polytope
P = struct; d = 10;
P.lb = zeros(d^2,1);
P.Aeq = sparse(2*d,d^2);
P.beq = ones(2*d,1);
for i=1:d
    P.Aeq(i,(i-1)*d+1:i*d) = 1;
    P.Aeq(d+i,i:d:d^2)=1;
end

opts = default_options();
opts.maxTime = 10; % Stop after 10 sec
fid = fopen('demo.log','w');
opts.outputFunc = @(tag, msg, varargin) {fprintf(fid, msg, varargin{:}); fprintf(msg, varargin{:})}; % output the results to the screen and demo.log
o = sample(P, +Inf, opts);
s = thin_samples(o.samples);
figure;
histogram(s(1,:))
title('Marginal of first coordinate of Birkhoff polytope');
fclose(fid);
drawnow()

%% Example 3: Sample Gaussian distribution restricted to a hypercube
P = struct; d = 100;
P.lb = -ones(d,1);
P.ub = ones(d,1);
P.f = @(x) x'*x/2;
P.df = @(x) x;
P.ddf = @(x) ones(d,1);
P.dddf = @(x) zeros(d,1);

opts = default_options();
opts.maxStep = 10000; % Stop after 10000 iter
o = sample(P, +Inf, opts);
s = thin_samples(o.samples);
figure;
histogram(s(:))
title('Marginal of Gaussian distribution restricted to hypercube');

%% Example 4: Read a polytope according to Cobra format
load('coverage\Recon1.mat')
P = struct;
P.lb = model.lb;
P.ub = model.ub;
P.S = model.b;
o = sample(P, 100);
[pVal] = uniformtest(o, struct('toPlot', true));

%% Example 5: Brownian bridge
P = struct; d = 10000;
e = ones(d,1);
P.Aeq = [spdiags([e -e], 0:1, d-1, d) spdiags(e, 0, d-1, d-1)];
P.beq = zeros(d-1,1);
P.lb = -10*ones(2*d-1,1);
P.ub = 10*ones(2*d-1,1);
P.lb(1:d) = -100 * sqrt(d);
P.ub(1:d) = 100 * sqrt(d);
P.lb([1 d]) = 0;
P.ub([1 d]) = 0;

P.f = @(x) x((d+1):end)'*x((d+1):end)/2;
P.df = @(x) [zeros(d,1);x((d+1):end)];
P.ddf = @(x) [zeros(d,1);ones(d-1,1)];
P.dddf = @(x) zeros(2*d-1,1);

o = sample(P, 10);
plot(o.samples(1:d,end))
title('Brownian bridge');
