#pragma once
#include <random>
#include "SparseMatrix.h"

// Problem:
// Approximate M = diag(A' inv(LL') A)

namespace PackedCSparse {
	template<typename Tx, typename Ti>
	struct LeverageJLOutput : DenseVector<Tx, Ti>
	{
		std::unique_ptr<Tx[]> d;		// random direction d
		std::unique_ptr<Tx[]> L_d;		// random direction d
		std::unique_ptr<Tx[]> AtL_d;	// A' L^{-1} d
		Ti k = 0, m = 0;
		std::mt19937_64 gen;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, Ti k_)
		{
			pcs_assert(L.initialized() && A.initialized() && At.initialized(), "leverageJL: bad inputs.");
			pcs_assert(k_ >= 1 && L.m == L.n && L.n == A.m && L.n == At.n && A.n == At.m, "leverageJL: dimensions mismatch.");

			if (k_ > k || A.m > this->m || A.n > this->n)
			{
				this->k = k_; this->n = A.n; this->m = A.m;
				this->x.reset(new Tx[this->n]);
				this->d.reset(new Tx[this->m * this->k]);
				this->L_d.reset(new Tx[this->m * this->k]);
				this->AtL_d.reset(new Tx[this->n * this->k]);
			}
		}
	};

	// compute 1/k sum_j (A' L^{-T} u_j) .* (A' L^{-T} v_j)
	// if diffProj, v and u are different random {+-1}. Otherwise v = u
	template<bool diffProj, size_t k, typename Tx, typename Ti, typename Tx2>
	void projectionJL(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		constexpr Ti k_ = k * (1 + diffProj);
		o.initialize(L, A, At, k_);

		Ti m = A.m, n = A.n;
		Tx T0 = Tx(0.0), T1 = Tx(1.0);
		Tx* d = o.d.get(), * L_d = o.L_d.get(), * AtL_d = o.AtL_d.get(), * x = o.x.get();

		for (Ti i = 0; i < m * k_; i++)
			d[i] = sign<Tx>(o.gen);

		for (Ti i = 0; i < n * k_; i++)
			AtL_d[i] = T0;

		ltsolve(L, (FloatArray<Tx, k_>*)d, (FloatArray<Tx, k_>*)L_d);
		gaxpy(At, (FloatArray<Tx, k_>*)L_d, (FloatArray<Tx, k_>*)AtL_d);

		Tx ratio = T1 / Tx(double(k));
		for (Ti i = 0; i < n; i++)
		{
			Tx ret_i = T0;
			for (Ti j = 0; j < k; j++)
				ret_i += AtL_d[i * k_ + j] * AtL_d[i * k_ + j + k * diffProj];

			x[i] = ratio * ret_i;
		}
	}

	template<size_t k, typename Tx, typename Ti, typename Tx2>
	void leverageJL(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		projectionJL<false, k>(o, L, A, At);
	}

	template<size_t k, typename Tx, typename Ti, typename Tx2>
	LeverageJLOutput<Tx, Ti> leverageJL(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		LeverageJLOutput<Tx, Ti> o;
		leverageJL<k>(o, L, A, At);
		return o;
	}

	template<size_t k, typename Tx, typename Ti, typename Tx2>
	Tx cholAccuracy(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, const Tx* w)
	{
		const Ti k_ = ((k & 1) == 1)?2*k:k; // I need even k so that FloatType<double, 2> is never called. 
		projectionJL<true, k_>(o, L, A, At);

		Ti m = A.m, n = A.n;
		Tx* x = o.x.get();

		Tx ret1 = Tx(0.0);
		for (Ti i = 0; i < n; i++)
			ret1 += x[i] * w[i];

		Tx ret2 = Tx(0.0);
		for (Ti i = 0; i < m; i++)
		{
			Tx* d = o.d.get() + i * (2 * k_);
			for (Ti j = 0; j < k_; j++)
				ret2 += d[j] * d[j + k_];

		}
		ret2 = ret2 / Tx(double(k_));

		Tx dist = (ret1 - ret2) * Tx(sqrt(double(k_)));
		return abs(dist);
	}
}