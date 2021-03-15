%% Example 5: Brownian bridge
n = 4;
dim = transpose(logspace(2, 5, n));
time = zeros(n, 1);
rh = zeros(n, 1);
ess = zeros(n, 1);
esspersec = zeros(n, 1);

for i = 1:n
    
    d = round(dim(i));
    P = struct; 
    e = ones(d,1);
    P.Aeq = [spdiags([e -e], 0:1, d-1, d) spdiags(e, 0, d-1, d-1)];
    P.beq = zeros(d-1,1);
    P.lb = -10*ones(2*d-1,1);
    P.ub = 10*ones(2*d-1,1);
    P.lb(1:d) = -100 * sqrt(d);
    P.ub(1:d) = 100 * sqrt(d);
    P.lb([1 d]) = 0;
    P.ub([1 d]) = 0;

    P.f = @(x) x((d+1):end)'*x((d+1):end)/2;
    P.df = @(x) [zeros(d,1);x((d+1):end)];
    P.ddf = @(x) [zeros(d,1);ones(d-1,1)];
    P.dddf = @(x) zeros(2*d-1,1);

    tic
    o = sample(P, 100);
    
    %s = thin_samples(o.samples(1:d, :));
    %plot(o.samples(1:d,end))
    %title('Brownian bridge');
    
    ess(i) = round(min(o.ess));
    time(i) = toc;
    rh(i) = max(abs(rhat(o.samples(1:d, :))-1)); 
    esspersec(i) = time(i)/ess(i);
end
esspersec = 1./esspersec;
table(dim, time, rh, ess, esspersec)