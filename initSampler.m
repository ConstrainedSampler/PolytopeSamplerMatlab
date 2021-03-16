function initSampler(recompile)
addpath(genpath(fullfile('code')));
addpath(fullfile('bin'));

addpath(fullfile('coverage'));
addpath(fullfile('coverage', 'problems'));

% check if the solver exists
if (~(nargin == 1 && recompile) && exist('PackedChol') == 3)
    try
        uid = PackedChol('init', uint64(123), speye(3));
        PackedChol('delete', uid);
        recompile = false;
    catch
        recompile = true;
    end
else
    recompile = true;
end

if (recompile)
    [path,~,~] = fileparts(mfilename('fullpath'));
    libpath = fullfile(path, 'code', 'solver');
    qdpath = fullfile(path, 'code', 'solver', 'qd');
    binpath = fullfile(path, 'bin');
    
    if (isempty(mex.getCompilerConfigurations('C++','Selected')))
        error('No C++ mex compiler is available.');
    end
    
    compiler = mex.getCompilerConfigurations('C++','Selected').ShortName;
    if (contains(compiler, 'MSVCPP'))
        mex_string = ...
            'mex -R2018a -silent -g -I"%libpath" COMPFLAGS="$COMPFLAGS /O2 /arch:AVX2"';
    elseif (contains(compiler, 'Clang++'))
        mex_string = ...
            'mex -R2018a -silent -O -I"%libpath" CFLAGS="$CFLAGS -march=native"';
    elseif (contains(compiler, 'g++'))
        mex_string = ...
            'mex -R2018a -silent -O -I"%libpath" CFLAGS="$CFLAGS -march=native"';
    else
        error('Currently, we only support MSVCPP, Clang++ or g++ as the compiler.');
    end
    
    c = {};
    c{end+1} = '%mex_string -output %binpath/PackedChol %libpath/PackedChol.cpp %qdpath/util.cc %qdpath/bits.cc %qdpath/dd_real.cc %qdpath/dd_const.cc %qdpath/qd_real.cc %qdpath/qd_const.cc';
    
    %% replace keywords
    keywords = {'%mex_string', '%qdpath', '%libpath', '%binpath'};
    replaces = {mex_string, qdpath, libpath, binpath};
    replaces = replace(replaces, keywords, replaces);

    for i = 1:length(c)
        c{i} = replace(c{i}, keywords, replaces);
    end

    %% Run the commands
    for i = 1:length(c)
        disp(c{i})
        eval(c{i});
    end
end