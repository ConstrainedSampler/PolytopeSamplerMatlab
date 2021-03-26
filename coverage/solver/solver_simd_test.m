function solver_simd_test(file, high_acc)
solver = @PackedChol4;

load(file);
A = problem.Aeq;
A = [A speye(size(A,1))];
w = rand(4, size(A,2)) + 0.2;
%w = ones(4,1) * rand(1, size(A,2)) + 0.2;
%w = 5 * ones(4, size(A,2));
uid = solver('init', uint64(1234), A);
if (high_acc)
    solver('setAccuracyTarget', uid, 0.0);
end
acc = solver('decompose', uid, w);

%% test L
L = solver('L', uid, 3);
H = A * diag(sparse(w(4,:))) * A';
assert(sum(abs(H - L * L'),'all') < 0.01)

%% test lsc
lsc = solver('leverageScoreComplement', uid, 0);

corank1 = sum(lsc, 2);
corank2 = size(A,2)-size(A,1);

assert(all(abs(corank1 - corank2)<0.01));
assert(all(lsc > -1e-5, 'all'));
assert(all(lsc < 1+1e-5, 'all'));

%% test lsc with JL
lsc_apx = solver('leverageScoreComplement', uid, 32);
assert(max(abs(lsc_apx-lsc), [], 'all') < 1.0);

%% test diagL
diagL = solver('diagL', uid);

for k = 1:4
L = chol(A * diag(sparse(w(k,:))) * A', 'lower');
diagL2 = diag(L);
assert(sum(abs(diagL(k,:)' - diagL2)) < 0.01);
end

%% test logdet
logdet = solver('logdet', uid);
logdet2 = sum(log(diagL), 2) * 2;

assert(all(abs(logdet - logdet2) < 0.01));

%% test solve
b = randn(4, size(A,1));
x = solver('solve', uid, b);
x2 = L'\ (L \ b(4,:)');

assert(sum(abs(x(4,:)' - x2)) < 0.01);

solver('delete', uid);
end