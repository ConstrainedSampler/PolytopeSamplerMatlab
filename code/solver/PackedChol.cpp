#include <cstring>
#include "mex_utils.h"
#include "PackedCSparse/PackedChol.h"

// SIMD_LEN = 0 means 1 and the vectors are inputted not as a matrix
#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif

namespace env = MexEnvironment;
typedef env::mexIdx Index;
const size_t chol_k = (SIMD_LEN == 0) ? 1 : SIMD_LEN;
using CholObj = PackedChol<chol_k, Index>;
using Matrix = SparseMatrix<double, Index>;
using Tx2 = FloatArray<double, chol_k>;

int main()
{
	size_t simd_len = SIMD_LEN;
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
		CholObj* solver = (CholObj*)uid;
		size_t n = solver->A.n; size_t m = solver->A.m;
		if (!strcmp(cmd, "solve"))
		{
			const double* b; double * out;
			if (SIMD_LEN == 0)
			{
				b = env::inputArray<double>(m);
				out = env::outputArray<double>(m);
			}
			else
			{
				b = env::inputArray<double>(simd_len, m);
				out = env::outputArray<double>(simd_len, m);
			}
			solver->solve((Tx2*)b, (Tx2*)out);
		}
		else if (!strcmp(cmd, "decompose"))
		{
			const double* w;
			if (SIMD_LEN == 0)
				w = env::inputArray<double>(n);
			else
				w = env::inputArray<double>(simd_len, n);
			Tx2 ret = solver->decompose((Tx2*)w);
			
			if (SIMD_LEN == 0)
				env::outputScalar<double>(get(ret, 0));
			else
			{
				double* out = env::outputArray<double>(simd_len, 1);
				for (size_t i = 0; i < SIMD_LEN; i++)
					out[i] = get(ret, i);
			}
		}
		else if (!strcmp(cmd, "leverageScoreComplement"))
		{
			size_t k = (size_t)env::inputScalar<double>();
			double* out;
			if (SIMD_LEN == 0)
				out = env::outputArray<double>(n);
			else
				out = env::outputArray<double>(simd_len, n);

			if (k == 0)
				solver->leverageScoreComplement((Tx2*)out);
			else
				solver->leverageScoreComplementJL((Tx2*)out, k);
		}
		else if (!strcmp(cmd, "logdet"))
		{
            Tx2 ret = solver->logdet();
			if (SIMD_LEN == 0)
				env::outputScalar<double>(get(ret, 0));
			else
			{
				double* out = env::outputArray<double>(simd_len, 1);
				for (size_t i = 0; i < SIMD_LEN; i++)
					out[i] = get(ret, i);
			}
		}
		else if (!strcmp(cmd, "diagL"))
		{
			double* out;
			if (SIMD_LEN == 0)
				out = env::outputArray<double>(m);
			else
				out = env::outputArray<double>(simd_len, m);
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
			throw "Invalid operation.";
	}
	mxFree((void*)cmd);
}