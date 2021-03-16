#include "mex_utils.h"
#include "PackedCSparse/PackedChol.h"

#define REPEAT_1 DECL(0)
#define REPEAT_2 DECL(0) DECL(1)
#define REPEAT_3 DECL(0) DECL(1) DECL(2)
#define REPEAT_4 DECL(0) DECL(1) DECL(2) DECL(3)
#define REPEAT_5 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4)
#define REPEAT_6 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5)
#define REPEAT_7 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6)
#define REPEAT_8 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7)
#define REPEAT_12 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11)
#define REPEAT_16 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15)
#define REPEAT_20 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19)
#define REPEAT_24 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23)
#define REPEAT_28 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23) DECL(24) DECL(25) DECL(26) DECL(27)
#define REPEAT_32 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23) DECL(24) DECL(25) DECL(26) DECL(27) DECL(28) DECL(29) DECL(30) DECL(31)

#define PROGRAM_REPEAT DECL(1) DECL(4) DECL(8) DECL(12) DECL(16) DECL(20) DECL(24) DECL(28) DECL(32)
#define SUPPORTED_MSG "We only support 1, 4, 8, 12, 16, 20, 24, 28 or 32 many vectors."


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
            size_t k = -1, m_ = m * chol_k;
            const double *b = env::inputArray<double>(k, m_);
			double *out = env::outputArray<double>(k, m_);
			
            #define DECL(k) case k: solver->solve<k>((Tx2*)b, (Tx2*)out); break;
            switch(k)
            {
                PROGRAM_REPEAT
                default:
                    env::error(SUPPORTED_MSG);
            }
            #undef DECL
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
            
            #define DECL(k) case k: solver->leverageScoreComplementJL<k>((Tx2*)out); break;
            switch(k)
            {
                PROGRAM_REPEAT
                case 0:
                    solver->leverageScoreComplement((Tx2*)out);
                    break;
                default:
                    env::error(SUPPORTED_MSG);
            }
            #undef DECL
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