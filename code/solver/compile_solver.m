function compile_solver(simd_len, recompile, arch)
if nargin < 1, simd_len = 0; end
if nargin < 2, recompile = false; end
if nargin < 3, arch = 'native'; end

% check if the solver exists
func_name = MexSolver.solverName(simd_len);
if (~recompile && exist(func_name) == 3)
   try
      pchol = str2func(func_name);
      uid = pchol('init', uint64(123), speye(3));
      pchol('delete', uid);
      recompile = false;
   catch
      recompile = true;
   end
else
   recompile = true;
end

if (recompile)
   clear mex
   [path,~,~] = fileparts(mfilename('fullpath'));
   libpath = fullfile(path);
   qdpath = fullfile(path, 'qd');
   outputfile = fullfile(path, '..', '..', 'bin', func_name);
   
   if (isempty(mex.getCompilerConfigurations('C++','Selected')))
      error('No C++ mex compiler is available.');
   end
   
   compiler = mex.getCompilerConfigurations('C++','Selected').ShortName;
   mex_string = ['mex -R2018a -silent -DSIMD_LEN=' num2str(simd_len) ' -O -I"%libpath" '];
   
   if (contains(func_name, 'arm'))
       archflag = '';
       nameflag = '';
       mex_string2 = ['CFLAGS="$CFLAGS ' archflag '"'];
   elseif (contains(compiler, 'MSVCPP'))
      switch arch
         case 'native'
            s = cpuInfo();
            if (s.AVX512F && s.OS_AVX512)
               archflag = '/arch:AVX512';
            elseif (s.AVX2 && s.OS_AVX)
               archflag = '/arch:AVX2';
            else
               archflag = '';
            end
            nameflag = 'native';
         case 'AVX'
            archflag = '/arch:AVX2';
            nameflag = '';
         case 'SSE'
            archflag = '';
            nameflag = 'SSE';
      end
      mex_string2 = ['COMPFLAGS="$COMPFLAGS /O2 ' archflag '"'];
   elseif (contains(compiler, 'Clang++'))
      switch arch
         case 'native'
            archflag = '-march=native';
            nameflag = 'native';
         case 'AVX'
            archflag = '-mavx2 -mfma';
            nameflag = '';
         case 'SSE'
            archflag = '';
            nameflag = 'SSE';
      end
      mex_string2 = ['CFLAGS="$CFLAGS ' archflag '"'];
   elseif (contains(compiler, 'g++'))
      switch arch
         case 'native'
            archflag = '-march=native';
            nameflag = 'native';
         case 'AVX'
            archflag = '-mavx2 -mfma';
            nameflag = '';
         case 'SSE'
            archflag = '';
            nameflag = 'SSE';
      end
      mex_string2 = ['CFLAGS="$CFLAGS ' archflag '"'];
   else
      error('Currently, we only support MSVCPP, Clang++ or g++ as the compiler.');
   end
   outputfile = [outputfile nameflag];
   mex_string = [mex_string mex_string2];
   
   c = {};
   c{end+1} = '%mex_string -output %outputfile %libpath/PackedChol.cpp %qdpath/util.cc %qdpath/bits.cc %qdpath/dd_real.cc %qdpath/dd_const.cc %qdpath/qd_real.cc %qdpath/qd_const.cc';
   
   %% replace keywords
   keywords = {'%mex_string', '%qdpath', '%libpath', '%outputfile'};
   replaces = {mex_string, qdpath, libpath, outputfile};
   replaces = replace(replaces, keywords, replaces);
   
   for i = 1:length(c)
      c{i} = replace(c{i}, keywords, replaces);
   end
   
   %% Run the commands
   for i = 1:length(c)
      disp(c{i});
      eval(c{i});
   end
end