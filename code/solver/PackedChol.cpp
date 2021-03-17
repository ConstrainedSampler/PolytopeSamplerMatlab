#include "mex_utils.h"
#include "PackedCSparse/PackedChol.h"

namespace env = MexEnvironment;
typedef env::mexIdx Index;
const size_t chol_k = 1;
using CholObj = PackedChol<chol_k, Index>;
using Matrix = SparseMatrix<double, Index>;
using Tx2 = FloatArray<double, chol_k>;

int main()
{
    const char* cmd = env::inputString();
    uint64_t uid = env::inputScalar<uint64_t>();
    if (!strcmp(cmd, "init"))
    {
        Matrix A = std::move(env::inputSparseArray<double>());
        auto solver = new CholObj(A);
        solver->setSeed(uid);
        env::outputScalar<uint64_t>((uint64_t)solver);
    }
    else
    {
        CholObj* solver = (CholObj *)uid;
        size_t n = solver->A.n; size_t m = solver->A.m;
        if (!strcmp(cmd, "solve"))
        {
            size_t m_ = m * chol_k;
            const double *b = env::inputArray<double>(m_);
			double *out = env::outputArray<double>(m_);
			solver->solve((Tx2*)b, (Tx2*)out);
        }
        else if (!strcmp(cmd, "decompose"))
        {
            size_t n_ = n * chol_k;
            auto w = env::inputArray<double>(n_);
            solver->decompose(w);
        }
        else if (!strcmp(cmd, "leverageScoreComplement"))
        {
            auto k = (size_t)env::inputScalar<double>(), n_ = n * chol_k;
            double *out = env::outputArray<double>(n_);
            
			if (k == 0)
				solver->leverageScoreComplement((Tx2*)out);
			else
				solver->leverageScoreComplementJL((Tx2*)out, k);
        }
        else if (!strcmp(cmd, "logdet"))
        {
            env::outputScalar<double>(solver->logdet());
        }
        else if (!strcmp(cmd, "diagL"))
        {
            auto out = env::outputArray<double>(m * chol_k);
            solver->diagL((Tx2*)out);
        }
        else if (!strcmp(cmd, "L"))
        {
            size_t k = (size_t)env::inputScalar<double>(0.0);
			auto L = solver->getL(k);
            env::outputSparseArray(L);
        }
        else if (!strcmp(cmd, "setAccuracyTarget"))
        {
            solver->accuracyThreshold = env::inputScalar<double>();
        }
        else if (!strcmp(cmd, "getDecomposeCount"))
        {
            env::outputDoubleArray(solver->numExact.data(), chol_k + 1);
        }
        else if (!strcmp(cmd, "delete"))
        {
            delete solver;
        }
        else
            env::error("Invalid operation.");
    }
}