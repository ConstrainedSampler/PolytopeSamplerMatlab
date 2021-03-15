debug = 1;
problems = [];
folder = {'basic', 'metabolic', 'netlib'};
%folder = {'basic', 'metabolic', 'netlib', 'extra'};

presolve_test(debug, folder, problems)
uniform_sample_test(debug, folder, problems)