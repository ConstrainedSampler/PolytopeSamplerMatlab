[path,~,~] = fileparts(mfilename('fullpath'));
file = fullfile(path, '..', 'problems', 'netlib', 'degen2');
load(file);
A = problem.Aeq;
A = [A 1e-8*speye(size(A,1))];
%rand('seed', 1);

for k = 1:1000
    rand('seed', k);
so1 = MultiMatlabSolver(A, 1e+100, 4);
so2 = MexSolver(A, 1e+100, 4);

w = rand(4, size(A,2)) + 0.2;
so1.setScale(w);
so2.setScale(w);

[so1.accuracy so2.accuracy]
%assert(norm(so1.diagL-so2.diagL) < 1)
assert(all([so1.accuracy so2.accuracy] < 1000000,'all'))
assert(~any(isnan([so1.accuracy so2.accuracy]),'all'))
end