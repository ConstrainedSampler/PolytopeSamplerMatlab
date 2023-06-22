debug = 0;
problems = [];
folders = {'basic', 'metabolic', 'netlib'};
presolve_test(debug, folders, problems)
%%

folders = {'basic', 'metabolic', 'netlib'};
sample_test(debug, folders, problems, 100, 'uniform')
sample_test(debug, folders, problems, 100, 'exponential')
sample_test(debug, folders, problems, 100, 'normal')


p_value_test();
