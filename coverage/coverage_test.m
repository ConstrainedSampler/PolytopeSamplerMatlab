debug = 0;
problems = [];
folders = {'basic', 'metabolic', 'netlib', 'extra', 'netlib_extra'};
presolve_test(debug, folders, problems)


folders = {'basic', 'metabolic', 'netlib'};
sample_test(debug, folders, problems, 1000, 'uniform')
sample_test(debug, folders, problems, 1000, 'expoential')
sample_test(debug, folders, problems, 1000, 'normal')


p_value_test();
