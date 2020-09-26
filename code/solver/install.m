c = {};

compiler = mex.getCompilerConfigurations('C++','Selected').ShortName;
if (contains(compiler, 'MSVCPP'))
    mex_with_fmath = ...
        'mex -R2018a -silent -g COMPFLAGS="$COMPFLAGS /O2 /arch:AVX2 /GL /fp:fast"';
    mex_no_fmath = ...
        'mex -R2018a -silent -g COMPFLAGS="$COMPFLAGS /O2 /arch:AVX2 /GL"';
    oext = '.obj';
elseif (contains(compiler, 'Clang++'))
    mex_with_fmath = ...
        'mex -R2018a -silent -O CFLAGS="$CFLAGS -march=native -ffast-math"';
    mex_no_fmath = ...
        'mex -R2018a -silent -O CFLAGS="$CFLAGS -march=native"';
    oext = '.o';
end

[path,~,~] = fileparts(mfilename('fullpath'));
obj = '%path/qd/util%oext %path/qd/bits%oext %path/qd/dd_real%oext %path/qd/dd_const%oext %path/qd/qd_real%oext %path/qd/qd_const%oext';

%% CSolver_double
c{end+1} = '%mex_fmath -output %path/CSolver_double %path/CSolver_double.cpp';

%% Compile qd objects
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/util.cc';
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/bits.cc';
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/dd_real.cc';
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/dd_const.cc';
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/qd_real.cc';
c{end+1} = '%mex_no_fmath -c -outdir %path/qd %path/qd/qd_const.cc';
c{end+1} = '%mex_no_fmath -output %path/CSolver_dd %path/CSolver_dd.cpp %obj';

%% replace keywords
keywords = {'%path', '%mex_fmath', '%mex_no_fmath', '%obj', '%oext'};
replaces = {path, mex_with_fmath, mex_no_fmath, obj, oext};
replaces = replace(replaces, keywords, replaces);

for i = 1:length(c)
    c{i} = replace(c{i}, keywords, replaces);
end

%% Run the commands
for i = 1:length(c)
    disp(c{i})
    eval(c{i});
end