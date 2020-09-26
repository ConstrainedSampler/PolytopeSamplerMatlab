function z = RungeKutta4(z, ham, opts)
% Perform the classical Runge–Kutta method
% This is 4-th order method, not symplectic.

%opts = setDefault(opts,default);

total_h = 0;
dz = ham.dH(z);
while total_h < 0.9*opts.trajLength
    ham.GenerateJL();
    h = min([   opts.maxRelativeStepSize * ham.StepSize(z, dz), ...
                opts.maxStepSize, ...
                opts.trajLength - total_h]);
    
    dz1 = ham.dH(z);
    
    h = min(h, ham.StepSize(z, dz1)/2);
    dz2 = ham.dH(z + h/2 * dz1);
    
    h = min(h, ham.StepSize(z, dz2)/2);
    dz3 = ham.dH(z + h/2 * dz2);
    
    h = min(h, ham.StepSize(z, dz3)/2);
    dz4 = ham.dH(z + h * dz3);
    dz = (dz1 + 2 * dz2 + 2 * dz3 + dz4) / 6;
    
    h = min(h, ham.StepSize(z, dz)/2);
    z = z + h * dz;
    
    if h < opts.minStepSize, z = NaN; return; end
    total_h = total_h + h;
end
end
