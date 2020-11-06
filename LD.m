function x = LD(N)
   
   x = zeros(1, N);
   eta = 0.01;
   for i = 2:N
       x(i) = x(i-1) - eta * 2 * x(i-1) + sqrt(2 * eta) * normrnd(0, 1);
   end
   
end