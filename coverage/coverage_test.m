debug = 0;
problems = [];
folders = {'basic', 'metabolic', 'netlib'};
%folders = {'basic', 'metabolic', 'netlib', 'extra'};

presolve_test(debug, folders, problems)
uniform_sample_test(debug, folders, problems, 200)