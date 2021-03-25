function update_script_ignore(debug)
if nargin == 0, debug = 0; end
initSampler_ignore

copyfile('C:/Users/yintat/OneDrive/MatrixLibrary/cpp/PackedCSparse', 'C:/GitHub/PolytopeSamplerMatlab/code/solver/PackedCSparse')
copyfile('C:/Users/yintat/OneDrive/MatrixLibrary/cpp/qd', 'C:/GitHub/PolytopeSamplerMatlab/code/solver/qd')
copyfile('C:/Users/yintat/OneDrive/MatrixLibrary/matlab/PackedChol.cpp', 'C:/GitHub/PolytopeSamplerMatlab/code/solver/PackedChol.cpp')
copyfile('C:/Users/yintat/OneDrive/MatrixLibrary/matlab/mex_utils.h', 'C:/GitHub/PolytopeSamplerMatlab/code/solver/mex_utils.h')

if (debug)
compile_solver_ignore(0, true);
compile_solver_ignore(1, true);
compile_solver_ignore(4, true);
else
compile_solver(0, true);
compile_solver(1, true);
compile_solver(4, true);
end

end