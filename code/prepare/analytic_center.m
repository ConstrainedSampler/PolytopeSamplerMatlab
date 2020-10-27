function [x, C, d, output] = analytic_center(A, b, f, opts, x)
%[x, C, d, output] = analytic_center(A, b, f, opts, x)
%compute the analytic center for the domain {Ax=b} intersect the domain of f
%
%Input:
% A - a m x n constraint matrix. The domain of the problem is given by Ax=b.
% b - a m x 1 constraint vector.
% f - a barrier class (Currently, we only have TwoSidedBarrier)
% opts - a structure for options with the following properties
%  display - whether we output the information at each iteration (default: false)
%  maxIter - maximum number of iterations (default: 300)
%  dualTol - stop when ||A' * lambda - gradf(x)||_2 < dualTol (default: 1e-12)
%  gaussianTerm - how much we add the identity into the hess f(x) to avoid numerical error (default: 1e-12)
%  regularizerStep - how much we multiply the identity we added into hess f(x) when there is a numerical problem (default: 10)
%  detectTightConstraints - detect tight constraints (default: true)
%  distanceTol - ||grad f_i(x)|| > 1/distanceTol && ||dx||_x >= velocityTol implies the coordinate i is tight. (default: 1e-6)
%  velocityTol - default: 1e-1.
%  distanceTol2 - ||grad f_i(x)|| > 1/distanceTol2 implies the coordinate i is tight. (default: 1e-9)
%
%Output:
% x - It outputs the minimizer of min f(x) subjects to {Ax=b}
% C - detected constraint matrix
%     If the domain ({Ax=b} intersect dom(f)) is not full dimensional in {Ax=b}
%     because of the dom(f), the algorithm will detect the collapsed dimension
%     and output the detected constraint C x = d
% d - detected constraint vector
% output - text output information at each iteration

%% prepare the printout
formats = struct;
formats.iter = struct('label', 'Iter', 'format', '5i');
formats.t = struct('label', 'Step Size', 'format', '13.2e');
formats.primalErr = struct('label', 'Primal Error', 'format', '13.2e');
formats.dualErr = struct('label', 'Dual Error', 'format', '13.2e');
output = TableDisplay(formats);
output.logFunc = opts.logFunc;
output.tag = 'analytic_center';
output.header();

%% initial conditions
if exist('x', 'var') == 0 || isempty(x) || ~f.barrier.feasible(x)
    x = f.barrier.center;
end
lambda = zeros(size(A,2),1);
fullStep = 0;
primalErr = Inf; dualErr = Inf; primalErrMin = Inf;
primalFactor = 1; dualFactor = 1 + norm(b);
idx = [];
solver = Solver(A, 'doubledouble');

%% find the central path
rx = lambda - f.gradient(x);
for iter = 1:opts.ipmMaxIter
    % compute the step direction
    rs = b - A * x;
    Hinv = 1./f.hessian(x);
    solver.setScale(Hinv);
    dr = solver.solve([rs A * (Hinv .* rx)]); % TODO
    dr1 = A' * dr(:,1); dr2 = A' * dr(:,2);
    dx1 = Hinv .* (dr1);
    dx2 = Hinv .* (rx - dr2);
    
    % compute the step size
    dx = dx1 + dx2;
    tGrad = min(f.barrier.step_size(x, dx),1);
    dx = dx1 + tGrad * dx2;
    tConst = min(0.99*f.barrier.step_size(x, dx),1);
    tGrad = tGrad * tConst;
    
    % update weight
    if opts.weightedBarrier
        if iter == 1
            f.barrier.weight = max(solver.leverageScoreComplement(),1e-6);
        else
            overage = solver.leverageScoreComplement() ./ f.barrier.weight;
            f.barrier.weight(overage > 2) = min(1, f.barrier.weight(overage > 2) * 2);
        end
    end
    
    if f.barrier.feasible(x + tConst * dx) % this is due to rounding error
        x = x + tConst * dx;
        lambda = lambda - dr2;
    
        % compute the residual
        rx = lambda - f.gradient(x);
        rs = b - A * x;
        primalErrMin = min(primalErr,primalErrMin); primalErr = norm(rx)/primalFactor;
        dualErrLast = dualErr; dualErr = norm(rs)/dualFactor;
    else
        dualErr = +Inf;
        primalErr = +Inf;
        dualErrLast = 0;
    end
    
    % check stagnation
    if ((dualErr > (1-0.9*tConst)*dualErrLast) && (primalErr > 10 * primalErrMin))
        % tight constraints condition: ||grad(x)||_2 > 1/distanceTol, (x is close boundary)
        dist = f.barrier.boundary_distance(x);
        idx = find(dist < opts.ipmDistanceTol);
    end
    
    % printout
    o = struct('iter', iter, 't', tGrad, 'primalErr', primalErr, 'dualErr', dualErr);
    output.print(o);
    
    % stopping criteria 
    if ~f.barrier.feasible(x) || ~isempty(idx)
        break;
    end
    
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
    note = sprintf('fixed %i barriers', length(idx));
    o = struct('iter', 0, 't', 0, 'primalErr', 0, 'dualErr', 0, 'note', note);
else
    C = zeros(0, size(A,2));
    d = zeros(0, 1);
end