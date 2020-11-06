P = struct; d = 10000;
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

o = sample(P, 100);
max(abs(rhat(o.samples(1:d, :))-1))
s = thin_samples(o.samples(1:d, :));
r = rhat(s);
max(abs(r-1))
plot(o.samples(1:d,end))
title('Brownian bridge');