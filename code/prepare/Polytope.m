classdef Polytope < handle
   % The problem is defined by exp(-f(z)) over {z = Tx+y where Ax = b, lb <= x <= ub}
   properties
      A       % constraint matrix A
      b       % constraint vector b
      T       % transformation T of the domain
      y       % shift of the domain
      barrier % TwoSidedBarrier
      f       % the objective function and its derivatives in the original space
      w      % Lewis weight
      df
      ddf
      original
      opts
      vdim	% oracles assumes each vector is along "dim"-th dimension
   end
   
   % derived variables
   properties
      center      % analytic center
      n
      width       % estimate of width of each variable
      
      fZero       % whether f is completely zero
      fHandle     % whether f is handle or not
      dfHandle    % whether df is handle or not
      ddfHandle   % whether ddf is handle or not
      
      % Assumed each row of T contains at most 1 non-zero
      Tidx        % T x = x(Tidx) .* Ta
      Ta          % T x = x(Tidx) .* Ta
      Tdf         % T' * df
      Tddf        % (T.^2)' * ddf
   end
   
   methods
      function o = Polytope(P, opts)
         % Input: a structure P with the following fields
         %  .Aineq
         %  .bineq
         %  .Aeq
         %  .beq
         %  .lb
         %  .ub
         %  .f
         %  .df
         %  .ddf
         % describing a log-concave distribution given by
         %   exp(-sum f_i(x_i))
         %       over
         %   {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}
         % where f is given by a vector function of its first and second
         % derivatives.
         %
         % Case 1: df is not defined
         %   f(x) = 0.
         % In this case, f, ddf must be empty.
         %
         % Case 2: df is a vector
         %   f_i(x_i) = df_i x_i.
         % In this case, f, ddf must be empty.
         %
         % Case 3: df is a function handle
         %   f need to be defined as a function handle.
         %   df need to be the derivative of f
         %   ddf is optional. Providing ddf improves the mixing time.
         if isa(P, 'Polytope')
            o = P;
            return;
         end
         P = standardize_problem(P);
         o.original = P;
         
         %% Convert the polytope into {Ax=b, lb<=x<=ub} form
         nP = size(P.Aeq, 2);
         nIneq = size(P.Aineq, 1);
         nEq = size(P.Aeq, 1);
         
         o.A = [P.Aeq sparse(nEq, nIneq); P.Aineq speye(nIneq)];
         o.b = [P.beq; P.bineq];
         o.f = P.f;
         o.df = P.df;
         o.ddf = P.ddf;
         lb = [P.lb; zeros(nIneq, 1)];
         ub = [P.ub; Inf*ones(nIneq, 1)];
         o.center = [];
         o.opts = opts.presolve;
         
         o.fHandle = isa(o.f, 'function_handle');
         o.dfHandle = isa(o.df, 'function_handle');
         o.ddfHandle = isa(o.ddf, 'function_handle');
         o.vdim = 1;
         
         if (~o.fHandle || ~o.dfHandle || ~o.ddfHandle)
            normf = norm(o.f) + norm(o.df) + norm(o.ddf);
            o.fZero = (normf == 0);
         end
         
         %% Move all variables with lb == ub to Ax = b
         fixedVars = dblcmp(ub, lb);
         if sum(fixedVars) ~= 0
            I = speye(length(lb));
            o.A = [I(fixedVars, :); o.A];
            o.b = [(lb(fixedVars)+ub(fixedVars))/2; o.b];
            ub(fixedVars) = +Inf;
            lb(fixedVars) = -Inf;
         end
         o.barrier = TwoSidedBarrier(lb, ub);
         o.barrier.extraHessian = opts.extraHessian;
         
         %% Update the transformation Tx + y
         o.T = sparse(nP, length(lb));
         o.T(:,1:nP) = speye(nP);
         o.updateT();
         o.y = zeros(nP, 1);
         
         %% Simplify the polytope
         if o.opts.runSimplify
            o.simplify();
            
            if isempty(o.center)
               o.opts.logFunc('Polytope:simplify', ['Run interior point methods to find the analytic center:' newline]);
               o.center = analytic_center(o.A, o.b, o, o.opts, o.center);
            end
            o.rescale(o.center);
            o.shift_barrier(o.center);
            o.reorder();
         else
            o.reorder();
         end
         
         %% Give Warning for Unbounded Polytope
         o.width = o.estimate_width(o.center);
         if (max(o.width) > 1e9)
            warning('Polytope:Unbounded', 'Domain seems to be unbounded. Either add a Gaussian term via f, df, ddf or add bounds to variable via lb and ub.');
         end
         
         %% Use the user-given center if it is given
         if ~isempty(P.center)
            o.center = o.T\(P.center - o.y);
         else
            %% Recenter again and make sure it is feasible
            [o.center, ~, ~, o.w] = lewis_center(o.A, o.b, o, o.opts, o.center);
            [~, hess] = o.lewis_center_oracle(o.center, o.w);
            solver = Solver(o.A, 'doubledouble');
            solver.setScale(1./hess);
            o.center = o.center + (o.A' * solver.solve(o.b - o.A*o.center))./hess;
            
            if (any(o.center > o.barrier.ub) || any(o.center < o.barrier.lb))
               error('Polytope:Infeasible', 'The algorithm cannot find a feasible point.');
            end
         end
      end
      
      function w = estimate_width(o, x)
         % Compute the width of the Dikin ellipse of the polytope:
         % w_i = sqrt(ei' H^{-1/2} (I - P) H^(-1/2} ei)
         
         solver = Solver(o.A, 'doubledouble');
         if nargin == 1 || isempty(x)
            hess = ones(size(o.A,2), 1);
         else
            [~, hess] = o.analytic_center_oracle(x);
         end
         solver.setScale(1./hess);
         tau = solver.leverageScoreComplement();
         w = sqrt((max(tau,0))./hess)+eps;
      end
      
      function simplify(o)
         o.opts.logFunc('Polytope:simplify', sprintf('Start simplifying the problem.\nInitially, there are %i variables and %i constraints.\n', o.n, size(o.A,1)));
         o.rescale();
         o.split_dense_cols(o.opts.splitDenseCols);
         o.reorder();
         o.remove_dependent_rows();
         changed = true;
         while changed
            while changed
               changed = false;
               changed = changed || o.remove_fixed_variables();
               changed = changed || o.remove_dependent_rows();
               o.reorder();
            end
            
            changed = o.extract_collapsed_variables();
         end
         o.opts.logFunc('Polytope:simplify', sprintf('Finish simplifying the problem.\nNow, there are %i variables and %i constraints.\n', o.n, size(o.A,1)));
      end
      
      function [g, h] = analytic_center_oracle(o, x)
         [~, g, h] = o.f_oracle(x);
         g = g + o.barrier.gradient(x);
         h = h + o.barrier.hessian(x);
      end
      
      function [g, h] = lewis_center_oracle(o, x, w)
         [~, g, h] = o.f_oracle(x);
         g = g + w .* o.barrier.gradient(x);
         h = h + w .* o.barrier.hessian(x);
      end
      
      function [f, g, h] = f_oracle(o, x)
         if (o.fZero)
            f = 0; g = 0; h = 0;
            return;
         end
         
         if o.vdim == 2
            x = x';
         end
         
         if (o.fHandle || o.dfHandle || o.ddfHandle)
            z = o.Ta .* x(o.Tidx,:) + o.y;
         end
         
         k = size(x, 2);
         if o.fHandle
            f = zeros(1, k);
            for i = 1:k
               f(i) = o.f(z(:,i));
            end
         else
            if (~o.dfHandle) % when f is just linear
               f = o.Tdf' * x;
            else
               f = o.f * ones(1,k);
            end
         end
         
         if o.dfHandle
            g = zeros(o.n, k);
            for i = 1:k
               g(o.Tidx, i) = o.Ta .* o.df(z(:,i));
            end
         else % when f is just linear
            g = o.Tdf * ones(1,k);
         end
         
         if o.ddfHandle
            h = zeros(o.n, k);
            for i = 1:k
               h(o.Tidx, i) = (o.Ta.*o.Ta) .* o.ddf(z(:,i));
            end
         else
            h = o.Tddf * ones(1,k);
         end
         
         if o.vdim == 2
            f = f';
            g = g';
            h = h';
         end
      end
      
      function set_vdim(o, vdim)
         o.vdim = vdim;
         o.barrier.set_vdim(vdim);
      end
   end
   
   methods (Access = private)
      function o = append_map(o, S, z)
         % Perform a change of variables
         % from the representation {Tx+y : x in P} to
         % {Tx+y, x = Sw+z : Sw+z in P}
         % (internal use only)
         
         % WARNING: This does not handle the ub and lb changes
         % It only update A, T and y
         if nargin < 3, z = zeros(size(S,1),1); end
         
         o.b = o.b - o.A * z;
         o.A = o.A * S;
         o.y = o.y + o.T * z;
         o.T = o.T * S;
         o.updateT();
      end
      
      function rescale(o, x)
         % Rescale the problem so it is in a better numerical form.
         
         % do not rescale if A has zero or one contraints/variables
         if min(size(o.A)) <= 1, return; end
         
         if nargin == 1 || isempty(x) || isempty(o.width)
            hess = ones(size(o.A,2), 1);
         else
            [~, hess] = o.analytic_center_oracle(x);
            hess = hess + 1./(o.width.*o.width);
         end
         
         [cscale,rscale] = gmscale(o.A./sqrt(hess)', 0, 0.9);
         o.A = spdiag(1./rscale) * o.A;
         o.b = o.b./rscale;
         o.barrier.set_bound(o.barrier.lb .* cscale, o.barrier.ub .* cscale);
         o.append_map(spdiag(1./cscale));
         
         if ~isempty(o.center)
            o.center = o.center .* cscale;
         end
      end
      
      function shift_barrier(o, x)
         o.append_map(speye(numel(x)), x);
         o.barrier.set_bound(o.barrier.lb - x, o.barrier.ub - x);
         
         if ~isempty(o.center)
            o.center = o.center - x;
         end
      end
      
      function changed = remove_fixed_variables(o)
         % Remove fixed variables implied by Ax = b
         tol = o.opts.removeFixedVariablesTol;
         
         % If the width < tol, we consider it as fixed variables
         d = o.estimate_width();
         solver = Solver(o.A, 'doubledouble');
         solver.setScale(ones(size(o.A,2),1));
         x = o.A'*solver.solve(o.b);
         x(abs(x) < tol) = 0;
         fixedVars = (d < tol*(1+abs(x)));
         if sum(fixedVars) > 0
            o.opts.logFunc('Polytope:simplify', sprintf('Removed %i fixed coordinates.\n', sum(fixedVars)));
            changed = true;
            % remove all fixed variables
            S = speye(o.n);
            o.append_map(S(:,~fixedVars), x.*fixedVars);
            o.barrier.set_bound(o.barrier.lb(~fixedVars), o.barrier.ub(~fixedVars));
            if ~isempty(o.center)
               o.center = o.center(~fixedVars);
            end
         else
            changed = false;
         end
      end
      
      function reorder(o)
         % Reorder vertices such that cholesky has better sparsity pattern
         
         % compute the cost of chol decomposition
         function s = cholCost(H, P)
            count = symbfact(H(P,P));
            s = sum(count.^2);
         end
         
         m = size(o.A,1);
         H = o.A * o.A' + spdiag(ones(m,1));
         
         p_dissect = dissect(H); p_amd = amd(H);
         
         if (cholCost(H, p_dissect) < cholCost(H, p_amd))
            p = p_dissect;
         else
            p = p_amd;
         end
         
         [~, Q] = etree(H(p,p));
         p = p(Q);
         
         o.A = o.A(p,:); o.b = o.b(p);
      end
      
      function changed = remove_dependent_rows(o)
         % Remove redundant rows in A
         
         % With zero row, diag(L) does not indicate dependent rows
         % correctly. So, we remove zero row first.
         zeroRows = full(sum(o.A~=0, 2))==0;
         o.A = o.A(~zeroRows,:); o.b = o.b(~zeroRows);
         
         solver = Solver(o.A, 'doubledouble');
         solver.setScale(ones(size(o.A,2),1))
         dL = solver.diagL();
         I = find((dL > 1e-12) .* (dL < 1e+64));
         
         if size(I,1) == size(o.A,1)
            changed = false;
            return; % This is important for empty matrix
         end
         
         o.opts.logFunc('Polytope:simplify', sprintf('Removed %i redundant constraints.\n', size(o.A,1) - size(I, 1)));
         o.A = o.A(I, :);
         o.b = o.b(I);
         changed = true;
      end
      
      function changed = extract_collapsed_variables(o)
         % Extract collapsed variables and move it to constraints
         o.opts.logFunc('Polytope:simplify', ['Run interior point methods to detect fixed coordinates:' newline]);
         [o.center, Ac, bc] = analytic_center(o.A, o.b, o, o.opts, o.center);
         
         changed = ~isempty(Ac);
         
         % update the A and b
         o.A = [Ac; o.A];
         o.b = [bc; o.b];
      end
      
      function split_dense_cols(o, maxNZ)
         % Rewrite P so that each cols has no more than maxNZ non-zeros
         
         if isempty(o.A) || size(o.A,1) <= maxNZ, return; end
         
         % return the original P if it is too dense
         if (nnz(o.A) > maxNZ * size(o.A,2))
            return;
         end
         
         A_ = o.A;
         nzCounts = full(sum(A_~=0,1));
         lb = o.barrier.lb; ub = o.barrier.ub;
         while max(nzCounts)>maxNZ
            [m_,n_] = size(A_);
            badCols = find(nzCounts > maxNZ);
            numBadCols = length(badCols);
            
            newA = cell(numBadCols, 1);
            for j = 1:numBadCols
               i = badCols(j);
               nzIndices = find(A_(:,i));
               midpoint = nzIndices(floor(length(nzIndices)/2));
               last = nzIndices(end);
               newA{j} = [zeros(midpoint,1);A_(midpoint+1:last,i);zeros(m_-last,1)];
               A_(midpoint+1:last,i) = 0;
            end
            
            A_ = horzcat(A_, newA{:});
            A_ = [A_; sparse(numBadCols, n_) -speye(numBadCols,numBadCols)];
            A_(m_+1:m_+numBadCols,badCols) = speye(numBadCols,numBadCols);
            o.b = [o.b; zeros(numBadCols,1)];
            lb = [lb; lb(badCols)];
            ub = [ub; ub(badCols)];
            
            nzCounts = full(sum(A_~=0,1));
         end
         
         o.T = [o.T sparse(size(o.T,1), length(ub)-size(o.T,2))];
         o.updateT();
         o.A = A_;
         o.n = length(lb);
         o.barrier.set_bound(lb, ub);
      end
      
      function updateT(o)
         % assume each row of T has at most 1 non-zeros
         assert(max(full(sum(o.T~=0,2))) <= 1);
         
         o.n = size(o.T,2);
         o.Tidx = int32((o.T~=0) * (1:o.n)');
         o.Tidx(o.Tidx==0) = 1;
         o.Ta = o.T * ones(o.n,1);
         
         if (~o.dfHandle) % when f is just linear
            if numel(o.df) == 0
               o.Tdf = 0;
            else
               o.Tdf = o.T' * o.df;
            end
         end
         
         if (~o.ddfHandle)
            if numel(o.ddf) == 0
               o.Tddf = 0;
            else
               o.Tddf = (o.T.*o.T)' * o.ddf;
            end
         end
      end
      
      
   end
end