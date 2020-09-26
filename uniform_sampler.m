% Sample uniform random points from a polytope. 
% Input: specification of polytope with equalities, inequalities and
% bounds on coordinates, N: number of samples.  
initSampler

%% Form the problem P
%name = 'metabolic/Recon3';
name = 'basic/random_sparse@10';
%name = 'netlib/etamacro';
P = loadProblem(name);
%special handling for netlib instances
if contains(name,'netlib/')
    x = linprog(P.df, P.Aineq, P.bineq, P.Aeq, P.beq, P.lb, P.ub, struct('Display','none'));
    threshold = P.df' * x + abs(P.df)' * abs(x);
    P.Aineq = [P.Aineq; P.df'];
    P.bineq = [P.bineq; threshold];
    P.ub = min(P.ub, 2 * max(abs(x)));
    P.lb = max(P.lb, -2 * max(abs(x)));
end
%uniform sampling in P
P.df = [];


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
