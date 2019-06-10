classdef Hamiltonian < handle
    properties
        % data
        A
        b
        f
        df
        barrier
        
        % JL options
        JLsize % 0 = not using JL.
        
        % internal variables
        JLdir
        m
        n
        so
    end
    methods
        function obj = Hamiltonian(barrier, A, b, f, df, JLsize)
            m = size(A,1); n = size(A,2);
            assert(all(size(b) == [m 1]));
            assert(barrier.n == n);
            
            obj.A = A;
            obj.b = b;
            obj.f = f;
            obj.df = df;
            obj.m = m;
            obj.n = n;
            obj.JLsize = JLsize;
            obj.barrier = barrier;
            obj.JLdir = [];
            obj.so = LinearSystemSolver(A);
            obj.GenerateJL();
        end
        
        function e = H(o,z)
            x = z(1:o.n);
            v = z((1+o.n):end);
            
            bar = o.barrier; A = o.A;
            
            hessInv = bar.HessianInv(x);
            o.so.Prepare(hessInv);
            g_inv_v = hessInv * v;
            IPv = g_inv_v - hessInv * (A' * o.so.Solve(A *g_inv_v));
            e = 0.5 * v' * IPv;
            
            R = o.so.getR();
            e = e + sum(log(diag(R)));
            H = bar.Hessian(x);
            e = e + sum(log(diag(H))) * 0.5;
            
            e = e + o.f(x);
        end
        
        function dz = dH(o, z)
            assert(all(size(z) == [2*o.n 1]));
            x = z(1:o.n);
            v = z((1+o.n):end);
            
            bar = o.barrier; A = o.A;
            if (any(isnan(z)) || ~bar.Feasible(x))
                dz = NaN;
                return;
            end
            
            hessInv = bar.HessianInv(x);
            o.so.Prepare(hessInv);
            g_inv_v = hessInv * v;
            
            zz = 6*(A*x-o.b) + A *g_inv_v;
            
            if (o.JLsize ~= 0)
                zz = [zz A *(bar.SqrtHessianInv(x)*o.JLdir)];
            end
            
            y = hessInv * (A' * o.so.Solve(zz));
                
            if (o.JLsize == 0)
                V = full((o.so.getR()'\(A * hessInv))');
            else
                %V = hessInv * (A' * o.so.Solve(A *(bar.SqrtHessianInv(x)*o.JLdir)));
                V = y(:,2:end);
            end
            sigma = bar.LogDetGradient(x) - bar.QuadraticFormGradient(x, V);
            
            %dx = g_inv_v - hessInv * (A' * o.so.Solve((A*x-o.b) + A *g_inv_v));
            dx = g_inv_v - y(:,1);
            dv = - o.df(x) + 0.5 * bar.QuadraticFormGradient(x, dx) - 0.5 * sigma;
            dz = [dx; dv];
        end
        
        function t = StepSize(o, z, dz)
            assert(all(size(z) == [2*o.n 1]));
            x = z(1:o.n);
            dx = dz(1:o.n);
            t1 = o.barrier.StepSize(x, dx);
            t2 = 1 / max(sqrt(o.barrier.HessianNorm(x, dx)));
            t = min(t1,t2);
        end
        
        function z = Generate(o, x)
            assert(all(size(x) == [o.n 1]));
            %o.JLdir = sign(randn(o.n,o.JLsize)) ./ sqrt(o.JLsize);
            
            % random direction
            gv = o.barrier.SqrtHessian(x) * randn(o.n,1);
            %o.so.Prepare(speye(o.n));
            %v = gv - o.A'*o.so.Solve(o.A*gv);
            
            % Above random direction has the same effect as the following:
            v = gv;
            %v = zeros(size(v));
            z = [x;v];
        end
        
        function GenerateJL(o)
            o.JLdir = sign(randn(o.n,o.JLsize)) ./ sqrt(o.JLsize);
        end
    end
end