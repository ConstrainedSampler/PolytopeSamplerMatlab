function t = timeit(f)
    t = tic;
    f(); % warm up
    t = toc(t);
    
    repeat = ceil(2 + 1e-5 / t); % the program takes at least 1e-5 sec.
    t = tic;
    for i = 1:repeat
        f();
    end
    t = toc(t) / repeat;
end
