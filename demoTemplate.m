%% Demo to sample
% Includes multiple ways to input instance, from a .mat file, from the
% problem library, define explicitly.
initSampler

%% Form the problem P
% Example 1: from a .mat file
% load('Recon2.v04.mat')
% P = struct;
% P.Aeq = modelR204.S;
% P.beq = modelR204.b;
% P.lb = modelR204.lb;
% P.ub = modelR204.ub;
% P.df = [];

% Examples 2,3: from the test library, uniform sampling
%P = loadProblem('basic/random_sparse@20');
%P.df = [];
%P = loadProblem('basic/birkhoff@100')
%P.df= [];

% Example 4: define directly (Gaussian in a cone)
% run scatter(out.samples(1,:),out.samples(2,:)) to see the projection to the first 2 coordinates
 P = struct;
 P.lb=zeros(1000,1);
 P.Aineq = zeros(1,size(P.lb,1));
 P.Aineq(1)=-1;
 P.Aineq(2)=-1;
 P.bineq = -1;
 P.df=@(x) x;
 P.f = @(x)x'*x/2;
 P.ddf = @(x) ones(size(P.lb,1),1);
 P.dddf = @(x) zeros(size(P.lb,1),1);

%% Samples from P

fid = fopen('log.txt','w');

profile on;
opts = default_options();
opts.maxTime = 3600*24;
%set to true for debugging info
opts.debugMode = false;
opts.nSketch = 0;
opts.minStepSize = 1e-10;
opts.dynamicBound = true;
opts.weightedBarrier = false;
opts.maxStep = 300000;
%if screen output not desired, use next line and comment out the following line. 
%opts.outputFunc = @(tag, msg, varargin) {fprintf(fid, msg, varargin{:})}; 
opts.outputFunc = @(tag, msg, varargin) {fprintf(fid, msg, varargin{:}); fprintf(msg, varargin{:})};

%Number of samples N
N=100;
o = sample(P, N, opts);

%profile report
if isempty(P.df)
    if opts.debugMode == true
        s = o.ham.T * o.samples + o.ham.y;
        [o.pVal] = uniformtest(s, P, o.dim,  struct('toPlot', 1));
    end
    [o.pVal] = uniformtest(o.samples, P, o.dim,  struct('toPlot', 1));
end
fclose(fid);
fprintf('See log.text for full output; samples in o.samples\n');
fprintf('Mixing Time: %fs,  Effective Sample Size: %i\n', o.mixingTime, round(min(o.ess)));

