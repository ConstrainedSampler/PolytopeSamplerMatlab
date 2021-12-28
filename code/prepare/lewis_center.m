function [x, C, d, wp] = lewis_center(A, b, f, opts, x)
%[x, C, d, wp] = lewis_center(A, b, f, opts, x)
%compute the lewis center for the domain {Ax=b} intersect the domain of f
%
%Input:
% A - a m x n constraint matrix. The domain of the problem is given by Ax=b.
% b - a m x 1 constraint vector.
% f - a barrier class (Currently, we only have TwoSidedBarrier)
% opts - a structure for options with the following properties
%  ipmMaxIter - maximum number of iterations
%  ipmDualTol - stop when ||A' * lambda - gradf(x)||_2 < dualTol
%  ipmDistanceTol - if x_i is ipmDistanceTol close to some boundary, we assume the coordinate i is tight.
%
%Output:
% x - It outputs the minimizer of min f(x) subjects to {Ax=b}
% C - detected constraint matrix
%     If the domain ({Ax=b} intersect dom(f)) is not full dimensional in {Ax=b}
%     because of the dom(f), the algorithm will detect the collapsed dimension
%     and output the detected constraint C x = d
% d - detected constraint vector

%% prepare the printout
formats = struct;
formats.iter = struct('label', 'Iter', 'format', '5i');
formats.t = struct('label', 'Step Size', 'format', '13.2e');
formats.primalErr = struct('label', 'Primal Error', 'format', '13.2e');
formats.dualErr = struct('label', 'Dual Error', 'format', '13.2e');
output = TableDisplay(formats);
opts.logFunc('analytic_center', output.header());

%% initial conditions
if exist('x', 'var') == 0 || isempty(x) || ~f.barrier.feasible(x)
    x = f.barrier.center;
end
lambda = zeros(size(A,2),1);
fullStep = 0; tConst = 0;
primalErr = Inf; dualErr = Inf; primalErrMin = Inf;
primalFactor = 1; dualFactor = 1 + norm(b);
idx = [];
solver = Solver(A, 'doubledouble');
w = ones(size(x));
wp = w;

%% find the central path
for iter = 1:opts.ipmMaxIter
    [grad, hess] = f.lewis_center_oracle(x, wp);
    
    % compute the residual
    rx = lambda - grad;
    rs = b - A * x;
    
    % check stagnation
    primalErrMin = min(primalErr,primalErrMin); primalErr = norm(rx)/primalFactor;
    dualErrLast = dualErr; dualErr = norm(rs)/dualFactor;
    feasible = f.barrier.feasible(x);
    if ((dualErr > (1-0.9*tConst)*dualErrLast) || (primalErr > 10 * primalErrMin) || ~feasible)
        dist = f.barrier.boundary_distance(x);
        idx = find(dist < opts.ipmDistanceTol);
        if ~isempty(idx), break; end
    end
    
    % compute the step direction
    Hinv = 1./hess;
    solver.setScale(Hinv);
    dr1 = A' * solver.solve(rs); dr2 = A' * solver.solve(A * (Hinv .* rx));
    dx1 = Hinv .* (dr1);
    dx2 = Hinv .* (rx - dr2);
    
    % compute the step size
    dx = dx1 + dx2;
    tGrad = min(f.barrier.step_size(x, dx),1);
    dx = dx1 + tGrad * dx2;
    tConst = min(0.99*f.barrier.step_size(x, dx),1);
    tGrad = tGrad * tConst;
    
    % make the step
    x = x + tConst * dx;
    lambda = lambda - dr2;
    
    % update weight
    wNew = max(double(solver.leverageScoreComplement(0)), 0) + 1e-8;%1 / length(w);
    w = (w + wNew)/2;
    wp = w.^(1-1/8);
    
    if ~f.barrier.feasible(x), break; end
    
    % printout
    o = struct('iter', iter, 't', tGrad, 'primalErr', primalErr, 'dualErr', dualErr);
    opts.logFunc('analytic_center', output.print(o));
    
    % stop if converged
    if (tGrad == 1)
        fullStep = fullStep + 1;
        if (fullStep > log(dualErr/opts.ipmDualTol) && fullStep > 8)
            break;
        end
    else
        fullStep = 0;
    end
end

if isempty(idx)
    dist = f.barrier.boundary_distance(x);
    idx = find(dist < opts.ipmDistanceTol);
end

if ~isempty(idx)
    [A_, b_] = f.barrier.boundary(x);
    C = A_(idx,:);
    d = b_(idx);
else
    C = zeros(0, size(A,2));
    d = zeros(0, 1);
end