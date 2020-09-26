function [x3, v3, step] = implicit_midpoint(x0, v0, h, ham, opts)
    % Step 1
    ham.prepare(x0);
    v1 = v0 - (h/2) * ham.DU(x0);
    x1 = ham.project(x0, h);
    done = 0;
    
    % Step 2
    x2 = x1; v2 = v1;
    nu = zeros(size(ham.A,1),1);
    for step = 1:opts.maxODEStep
        x2_old = x2;
        xmid = (x1+x2)/2;
        vmid = (v1+v2)/2;
        
        [dKdv, dKdx, nu] = ham.approxDK(xmid, vmid, nu);
        x2 = x1 + h * dKdv;
        v2 = v1 - h * dKdx;
        
        dist = ham.x_norm(xmid, x2-x2_old)/ h;
        if (dist < opts.implicitTol)
            done = 1;
            break;
        elseif any(~isfinite(dist))% || dist > 1e8
            break;
        end
    end
    
    if done == 0
        x3 = NaN; v3 = NaN;
        return
    end
    
    % Step 3
    x3 = x2;
    v3 = v2 - (h/2) * ham.DU(x3);
end
