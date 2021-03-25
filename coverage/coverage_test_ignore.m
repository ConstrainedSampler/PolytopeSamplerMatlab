debug = 1;
problems = [];
%folders = {'metabolic'};
%folders = {'basic', 'metabolic', 'netlib'};
folders = {'netlib'};
%folders = {'basic', 'metabolic', 'netlib', 'extra'};

%presolve_test(debug, folders, problems)
uniform_sample_test(debug, folders, problems, 20)
%p_value_test();