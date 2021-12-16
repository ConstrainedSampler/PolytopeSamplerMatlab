function RHMC_compile_all(internal_use)
[path, ~, ~] = fileparts(mfilename('fullpath'));
outputfile = fullfile(path, '..', '..', 'bin', 'cpuInfo');
mex_string = ['mex ' fullfile(path, 'FeatureDetector', 'cpuInfo.cpp') ' -output ' outputfile];
disp(mex_string);
eval(mex_string);

%%
if (nargin == 1 && internal_use)
   arch = 'AVX';
   compile_solver(0, true, arch);
   compile_solver(1, true, arch);
   compile_solver(4, true, arch);
   arch = 'SSE';
   compile_solver(0, true, arch);
   compile_solver(1, true, arch);
   compile_solver(4, true, arch);
else
   arch = 'native';
   compile_solver(0, true, arch);
   compile_solver(1, true, arch);
   compile_solver(4, true, arch);
end
