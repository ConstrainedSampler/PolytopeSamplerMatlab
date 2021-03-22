function solver_test()
[path,~,~] = fileparts(mfilename('fullpath'));
matrix_file = fullfile(path, '..', 'problems', 'netlib', 'degen2');
solver_scalar_test(matrix_file, false);
solver_scalar_test(matrix_file, true);
solver_zero_test(true, 0);
solver_zero_test(true, 4);
solver_zero_test(false, 0);
solver_zero_test(false, 4);
solver_simd_test(matrix_file, false);
solver_simd_test(matrix_file, true);
