% Solve the equation o.AwAt y = b with initial point x
% with an approximate solver o.approxSolve
%
% flag = 0 means success
% flag = 1 means fails
% flag = 2 means haven't converged
function [x, iter, flag] = batch_pcg(x, b, o, tol, max_iter)
    if (o.k == 0)
        vdim = 1;
    else
        vdim = 2;
    end
    
    flag = 0; % success
    r = b - o.AwAt(x);
    rz = 1;
    b_norm = sqrt(sum(b.^2,vdim));
    r_norm = sqrt(sum(r.^2,vdim));

    if all(r_norm <= tol * b_norm)
        iter = 0;
        return
    end

    for iter = 1:max_iter
        z = o.approxSolve(r);

        rz_prev = rz;
        rz = sum(r.*z,vdim);

        if any(rz < -eps) || any(~isfinite(rz))
            flag = 1; % stagnation
            break;
        end

        if iter == 1
            p = z;
        else
            beta = rz ./ rz_prev;
            p = z + beta .* p;
        end

        Ap = o.AwAt(p);
        pAp = sum(p.*Ap,vdim);
        if any(pAp <= 0) || any(~isfinite(pAp))
            flag = 1; % stagnation
            break;
        end
        alpha = rz./sum(p.*Ap,vdim);

        x = x + alpha .* p;
        r = r - alpha .* Ap;

        r_norm = sqrt(sum(r.^2,vdim));

        if all(r_norm <= tol * b_norm)
            return;
        end
    end
    
    if (flag ~= 1), flag = 2; end
    x = o.approxSolve(b);
end 