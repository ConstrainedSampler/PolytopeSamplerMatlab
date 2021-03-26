function solver_scalar_test(file, high_acc)
solver = @PackedChol0;

load(file);
A = problem.Aeq;
A = [A speye(size(A,1))];
w = ones(size(A,2),1);
uid = solver('init', uint64(1234), A);
if (high_acc)
    solver('setAccuracyTarget', uid, 0.0);
end
acc = solver('decompose', uid, w);

%% test lsc
lsc = solver('leverageScoreComplement', uid, 0);

corank1 = sum(lsc);
corank2 = size(A,2)-size(A,1);

assert(abs(corank1 - corank2)<0.01);
assert(all(lsc > -1e-5));
assert(all(lsc < 1+1e-5));

%% test lsc with JL
lsc_apx = solver('leverageScoreComplement', uid, 32);
assert(max(abs(lsc_apx-lsc)) < 1.0);

%% test diagL
diagL = solver('diagL', uid);

L = chol(A * diag(sparse(w)) * A', 'lower');
diagL2 = diag(L);
assert(sum(abs(diagL - diagL2)) < 0.01);

%% test logdet
logdet = solver('logdet', uid);
logdet2 = sum(log(diagL2)) * 2;

assert(abs(logdet - logdet2) < 0.01);

%% test solve
b = randn(size(A,1),1);
x = solver('solve', uid, b);
x2 = L'\ (L \ b);

assert(sum(abs(x - x2)) < 0.01);


solver('delete', uid);
end