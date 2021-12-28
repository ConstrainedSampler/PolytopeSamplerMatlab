function [x4, v4, step] = gauss_legendre(x0, v0, h, ham, opts)
    % Step 1
    x1 = x0;
    v1 = v0 - (h/2) * ham.DU(x0);
    done = 0;

    % Step 2
    c1 = h/4; c2 = h * (1/4-sqrt(3)/6);
    c3 = h * (1/4+sqrt(3)/6); c4 = h/4;
    nu1 = zeros(size(x1,1),size(ham.A,1),1); nu2 = zeros(size(x1,1),size(ham.A,1),1);
    k1x = zeros(size(x1)); k1v = zeros(size(v1));
    k2x = zeros(size(x1)); k2v = zeros(size(v1));
    
    for step = 1:opts.maxODEStep
        k2x_old = k2x;
        
        [k1x, k1v, nu1] = ham.approxDK...
            (x1 + c1 * k1x + c2 * k2x, ...
            v1 - c1 * k1v - c2 * k2v, nu1);
        
        [k2x, k2v, nu2] = ham.approxDK...
            (x1 + c3 * k1x + c4 * k2x, ...
            v1 - c3 * k1v - c4 * k2v, nu2);
        
        dist = ham.x_norm(x1, k2x_old-k2x);
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
    
    x2 = x1 + h/2 * (k1x + k2x);
    v2 = v1 - h/2 * (k1v + k2v);
    
    % Step 3
    x3 = x2;
    v3 = v2 - (h/2) * ham.DU(x3);
    
    % Step 4 (Project to Ax = b)
    ham.prepare(x3);
    v4 = v3;
    x4 = ham.project(x3);
end
