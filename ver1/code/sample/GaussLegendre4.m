function z = GaussLegendre4(z, ham, opts)
% Perform the collocation method based on Gauss Legendre with 2 points
% This is 4-th order symplectic method
% 
% options:
%  implicitIter - number of fix point iteration in the implicit step (default:2)

default.implicitIter = 2;
opts = setDefault(opts,default);

total_h = 0;
dz = ham.dH(z);
while total_h < 0.9*opts.trajLength
    ham.GenerateJL();
    h = min([   opts.maxRelativeStepSize * ham.StepSize(z, dz), ...
                opts.maxStepSize, ...
                opts.trajLength - total_h]);
    
    dz1 = dz; dz2 = dz;
    for i = 1:opts.implicitIter
        dz1_ = 0.25 * dz1 + (1/4-sqrt(3)/6) * dz2;
        h = min(h, ham.StepSize(z, dz1_)/2);
        dz1 = ham.dH(z + h * dz1_);
        
        dz2_ = (1/4+sqrt(3)/6) * dz1 + 1/4 * dz2;
        h = min(h, ham.StepSize(z, dz2_)/2);
        dz2 = ham.dH(z + h * dz2_);
    end
    
    dz = (dz1+dz2)/2;
    z = z + h * dz;
    
    if h < opts.minStepSize, z = NaN; return; end
    total_h = total_h + h;
end
end