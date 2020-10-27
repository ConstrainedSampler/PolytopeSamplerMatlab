classdef Polytope < handle
    % The problem is defined by exp(-f(x)) over {Tx+y where Ax=b, lb <= x <= ub}
    properties
        A   % constraint matrix A
        b   % constraint vector b
        T	% transformation T of the domain
        y	% shift of the domain
        barrier % TwoSidedBarrier
        f	% the objective function and its derivatives in the original space
        df
        ddf
        %dddf
        originalProblem
        
        opts
        
        % private variables
        center      % analytic center
        n
        
        % TODO: Assumed each row of T contains at most 1 non-zero
        T2          % used for computing ddf
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
            %  .dddf
            % describing a log-concave distribution given by
            %   exp(-sum f_i(x_i))
            %       over
            %   {Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}
            % where f is given by a vector function of its 1-st, 2-nd, 3-rd
            % derivative.
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
            %   df need to be the deriative of f
            %   ddf is optional. Providing ddf improves the mixing time.
            if isa(P, 'Polytope')
                o = P;
                return;
            end
            P = standardize_problem(P);
            o.originalProblem = P;
            
            %% Convert the polytope into {Ax=b, lb<=x<=ub} form
            nP = size(P.Aeq, 2);
            nIneq = size(P.Aineq, 1);
            nEq = size(P.Aeq, 1);
            
            o.A = [P.Aeq sparse(nEq, nIneq); P.Aineq speye(nIneq)];
            o.b = [P.beq; P.bineq];
            o.f = P.f;
            o.df = P.df;
            o.ddf = P.ddf;
            %o.dddf = P.dddf;
            lb = [P.lb; zeros(nIneq, 1)];
            ub = [P.ub; Inf*ones(nIneq, 1)];
            o.center = [];
            o.n = length(lb);
            o.opts = opts;
            
            %% Move all variables with lb == ub to Ax = b
            fixedVars = dblcmp(ub, lb);
            if sum(fixedVars) ~= 0
                I = speye(o.n);
                o.A = [I(fixedVars, :); o.A];
                o.b = [(lb(fixedVars)+ub(fixedVars))/2; o.b];
                ub(fixedVars) = +Inf;
                lb(fixedVars) = -Inf;
                o.n = length(lb);
            end
            o.barrier = TwoSidedBarrier(lb, ub);
            o.barrier.extraHessian = o.opts.extraHessian;
            
            %% Update the transformation Tx + y
            o.T = sparse(nP, o.n);
            o.T(:,1:nP) = speye(nP);
            o.T2 = o.T.^2;
            o.y = zeros(nP, 1);
            
            %% Simplify the polytope
            if o.opts.runSimplify
                o.simplify();
            else
                o.reorder();
            end
            
            %% Find an initial point
            if ~isempty(P.center)
                o.center = P.center;
            elseif ~isempty(o.center)
                o.opts.weightedBarrier = false;
                o.center = analytic_center(o.A, o.b, o, o.opts, o.barrier.center);
            end
            
            %% Give Warning for Unbounded Polytope
            w = o.estimate_width(o.center);
            if (max(w) > 1e10)
                warning('Hamiltonian:Unbounded', 'Domain seems to be unbounded. Either add a Gaussian term via f, df, ddf or add bounds to variable via lb and ub.');
            end
            
            %% Compute the initial weight for weighted barrier
%             if opts.weightedBarrier
%                 o.barrier = WeightedTwoSidedBarrier(o.barrier.lb, o.barrier.ub, ones(o.n,1));
%                 o.barrier.extraHessian = o.opts.extraHessian;
%                 o.opts.weightedBarrier = true;
%                 o.center = analytic_center(o.A, o.b, o, o.opts, o.center);
%             end
        end
        
        function w = estimate_width(o, x)
            % Compute the width of the Dikin ellipse of the polytope:
            % w_i = sqrt(ei' H^{-1/2} (I - P) H^(-1/2} ei)
            
            solver = Solver(o.A, 'doubledouble');
            if nargin == 1 || isempty(x)
                hess = ones(size(o.A,2), 1);
            else
                hess = o.hessian(x);
            end
            solver.setScale(1./hess);
            tau = solver.leverageScoreComplement();
            w = sqrt(tau./hess)+eps;
        end
        
        function simplify(o)
            o.rescale();
            o.split_dense_cols(o.opts.splitDenseCols);
            o.reorder();
            nTmp = +Inf;
            o.remove_dependent_rows();
            o.remove_fixed_variables();
            
            while o.n < nTmp
                nTmp = o.n;
                o.extract_collapsed_variables();
                o.remove_dependent_rows();
                o.remove_fixed_variables();
                
            end
            o.reorder();
        end
        
        function grad = gradient(o, x)
            grad = o.barrier.gradient(x) + o.T' * o.df(o.T * x + o.y);
        end
        
        function h = hessian(o, x)
            h = o.barrier.hessian(x) + o.T2' * o.ddf(o.T * x + o.y);
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
            o.T2 = o.T.^2;
        end
        
        function rescale(o)
            % Rescale the problem so it is in a better numerical form.
            
            % do not rescale if A has zero or one contraints/variables
            if min(size(o.A)) <= 1, return; end
            
            [cscale,rscale] = gmscale(o.A, 0, 0.9);
            o.A = spdiag(1./rscale) * o.A;
            o.b = o.b./rscale;
            o.barrier.update(o.barrier.lb .* cscale, o.barrier.ub .* cscale);
            o.append_map(spdiag(1./cscale));
        end
        
        function remove_fixed_variables(o)
            % Remove fixed variables implied by Ax = b
            tol = o.opts.removeFixedVariablesTol;
            
            % If the width < tol, we consider it as fixed variables
            d = o.estimate_width();
            solver = Solver(o.A, 'doubledouble');
            solver.setScale(ones(size(o.A,2),1));
            x = o.A'*solver.solve(o.b);
            x(abs(x) < 1e-8) = 0; 
            fixedVars = (d < tol*(1+abs(x)));
            
            % remove all fixed variables
            S = speye(o.n);
            o.append_map(S(:,~fixedVars), x.*fixedVars);
            o.barrier.update(o.barrier.lb(~fixedVars), o.barrier.ub(~fixedVars));
            if ~isempty(o.center)
                o.center = o.center(~fixedVars);
            end
            o.n = length(o.barrier.lb);
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
        
        function remove_dependent_rows(o)
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
                return; % This is important for empty matrix
            end
            
            o.A = o.A(I, :);
            o.b = o.b(I);
        end
        
        function extract_collapsed_variables(o)
            % Extract collapsed variables and move it to constraints
            o.opts.weightedBarrier = false;
            [o.center, Ac, bc] = analytic_center(o.A, o.b, o, o.opts, o.center);
            
            % update the A and b
            o.A = [Ac; o.A];
            o.b = [bc; o.b];
        end
        
        function split_dense_cols(o, maxNZ)
            % Rewrite P so that each cols has no more than maxNZ non-zeros
            
            if isempty(o.A) || size(o.A,1) <= maxNZ, return; end
            
            % return the original P if it is too dense
            if (nnz(o.A) > maxNZ * size(o.A,1))
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
            o.T2 = o.T.^2;
            o.A = A_;
            o.n = length(lb);
            o.barrier.update(lb, ub);
        end
    end
end