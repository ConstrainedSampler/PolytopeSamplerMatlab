%% Example 1: Sample uniform from a simplex
initSampler

P = struct; d = 100;
P.Aineq = ones(1, d);
P.bineq = 1;
P.lb = zeros(d, 1);
o = sample(P, 2000); % Number of samples = 2000
s = o.samples;
histogram(sum(s), 0.9:0.005:1)
title('distribution of l1 norm of simplex');
uniformtest(o, struct('toPlot', true));
drawnow()

%% Example 2: Sample uniform from Birkhoff polytope
initSampler

P = struct; d = 10;
P.lb = zeros(d^2,1);
P.Aeq = sparse(2*d,d^2);
P.beq = ones(2*d,1);
for i=1:d
    P.Aeq(i,(i-1)*d+1:i*d) = 1;
    P.Aeq(d+i,i:d:d^2)=1;
end

opts = default_options();
opts.maxTime = 20; % Stop in 20 sec
opts.logging = 'demo_ignore.log'; % Output the debug log to demo_ignore.log
o = sample(P, +Inf, opts);
s = o.samples;
figure;
histogram(s(1,:), 'BinLimits', [0, 0.5])
title('Marginal of first coordinate of Birkhoff polytope');
drawnow()

%% Example 3: Sample Gaussian distribution restricted to a hypercube
initSampler

P = struct; d = 1000;
P.lb = -ones(d,1);
P.ub = ones(d,1);
P.f = @(x) x'*x/2;
P.df = @(x) x;
P.ddf = @(x) ones(d,1);

opts = default_options();
opts.maxStep = 10000; % Stop after 10000 iter
o = sample(P, +Inf, opts);
figure;
histogram(o.samples(:), 'BinLimits', [-1, 1])
title('Marginal of Gaussian distribution restricted to hypercube');
drawnow()

%% Example 4: Brownian bridge
initSampler

P = struct; d = 1000;
e = ones(d,1);
P.Aeq = [spdiags([e -e], 0:1, d-1, d) spdiags(e, 0, d-1, d-1)];
P.beq = zeros(d-1,1);
P.lb = -Inf*ones(2*d-1,1);
P.ub = Inf*ones(2*d-1,1);
P.lb([1 d]) = 0;
P.ub([1 d]) = 0;

P.f = @(x) x((d+1):end)'*x((d+1):end)/2;
P.df = @(x) [zeros(d,1);x((d+1):end)];
P.ddf = @(x) [zeros(d,1);ones(d-1,1)];

o = sample(P, 100);
figure;
plot(o.samples(1:d, end))
title('Brownian bridge');

%% Example 5: Read a polytope according to Cobra format
initSampler

load(fullfile('coverage','Recon1.mat'))
P = struct; % Warning: Other Cobra models may have optional constraints (C,d)
P.lb = model.lb;
P.ub = model.ub;
P.beq = model.b;
P.Aeq = model.S;
o = sample(P, 100);
uniformtest(o, struct('toPlot', true));

%% Example 6: Run the sampler in parallel
initSampler
if isempty(ver('parallel'))
    fprintf('Parallel Computing Toolbox is required for this example')
else
    load(fullfile('coverage','Recon1.mat'))
    P = struct;
    P.lb = model.lb;
    P.ub = model.ub;
    P.beq = model.b;
    P.Aeq = model.S;
    opts = default_options();
    opts.nWorkers = 0;  % 0 means the default number of workers in the Parallel Computing Toolbox
    o = sample(P, 200, opts);
    [pVal] = uniformtest(o, struct('toPlot', true));
end