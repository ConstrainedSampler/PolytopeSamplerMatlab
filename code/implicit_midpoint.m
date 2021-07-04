function [x4, v4, step] = implicit_midpoint(x0, v0, h, ham, opts)
    % Step 1
    x1 = x0;
    v1 = v0 - (h/2) * ham.DU(x0);
    done = 0;
    
    % Step 2
    x2 = x1; v2 = v1;
    nu = zeros(size(x1,1),size(ham.A,1),1);
    for step = 1:opts.maxODEStep
        x2_old = x2;
        xmid = (x1+x2)/2;
        vmid = (v1+v2)/2;
        [dKdv, dKdx, nu] = ham.approxDK(xmid, vmid, nu);
        x2 = x1 + h * dKdv;
        v2 = v1 - h * dKdx;
        dist = ham.x_norm(xmid, x2-x2_old)/ h;
        if (max(dist,[],'all') < opts.implicitTol)
            done = 1;
            break;
        elseif any(dist > 1e16, 'all')
            break;
        end
    end
    
    if done == 0
        x4 = NaN; v4 = NaN;
        return
    end
    
    % Step 3
    x3 = x2;
    v3 = v2 - (h/2) * ham.DU(x3);
    
    % Step 4 (Project to Ax = b)
    ham.prepare(x3);
    v4 = v3;
    x4 = ham.project(x3);
end
