function solver_zero_test(high_acc, simd_len)
solver = str2func(['PackedChol' num2str(simd_len)]);

A = sparse(5,5);

if simd_len == 0
    w = ones(size(A,2),1);
else
    w = ones(simd_len, size(A,2));
end

uid = solver('init', uint64(1234), A);
if (high_acc)
    solver('setAccuracyTarget', uid, 0.0);
end
acc = solver('decompose', uid, w);
x = solver('leverageScoreComplement', uid, 0);
x = solver('leverageScoreComplement', uid, 32);
x = solver('diagL', uid);
x = solver('L', uid);
x = solver('logdet', uid);
if simd_len == 0
    b = randn(size(A,1),1);
else
    b = randn(simd_len, size(A,1));
end
x = solver('solve', uid, b);