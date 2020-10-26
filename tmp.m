x = rand(10000,1);
y = rand(10000,1);
z = x.*y;
t = z - z .* log(z);
hist(t)

x = 0:0.0001:1;
t = -x./lambertw(-1,-x/exp(1));
plot(x,t);


x = 0:0.0001:1;
t = (x) .* log(x/exp(1));

plot(x,t);