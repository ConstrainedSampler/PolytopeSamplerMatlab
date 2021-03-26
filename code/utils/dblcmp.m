% compares a and b and returns true if they are identical up to tol
function c = dblcmp(a, b, tol)
    if nargin < 3, tol = eps; end
    
    % <= does not work for the case with Inf
    c = (abs(a - b) < tol * (1+abs(a)+abs(b)));
end