nWorkers = 6;
p = gcp();
if nWorkers ~= 0 && p.NumWorkers ~= nWorkers
    delete(p);
    p = parpool(nWorkers);
end
nWorkers = p.NumWorkers;

x = 5;
spmd
    A = randn(1000);
    for i = 1:100
        if labindex == 1
            disp(sprintf('\b\b\b\b\b\b\b\b\b%6i', i));
        end
        B = svd(A);
    end
end