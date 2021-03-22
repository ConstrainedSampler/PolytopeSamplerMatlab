#pragma once
#include <random>
#include "SparseMatrix.h"

// Problem:
// Approximate M = diag(A' inv(LL') A)
namespace PackedCSparse {
	const size_t JLPackedSize = 4;

	template<typename Tx, typename Ti>
	struct LeverageJLOutput : DenseVector<Tx, Ti>
	{
		UniqueAlignedPtr<Tx> d;		// random direction d
		UniqueAlignedPtr<Tx> L_d;		// random direction d
		UniqueAlignedPtr<Tx> AtL_d;	// A' L^{-1} d
		Ti m = 0;
		std::mt19937_64 gen;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
		{
			pcs_assert(L.initialized() && A.initialized() && At.initialized(), "leverageJL: bad inputs.");
			pcs_assert(L.m == L.n && L.n == A.m && L.n == At.n && A.n == At.m, "leverageJL: dimensions mismatch.");
			this->n = A.n; this->m = A.m;
			this->x.reset(pcs_aligned_new<Tx>(this->n));
			this->d.reset(pcs_aligned_new<Tx>(this->m * 2 * JLPackedSize));
			this->L_d.reset(pcs_aligned_new<Tx>(this->m * 2 * JLPackedSize));
			this->AtL_d.reset(pcs_aligned_new<Tx>(this->n * 2 * 2 * JLPackedSize));
		}
	};

	// compute sum_{j=1}^{k} (A' L^{-T} u_j) .* (A' L^{-T} v_j)
	// if diffProj, v and u are different random {+-1}. Otherwise v = u
	template<bool diffProj, size_t k, typename Tx, typename Ti, typename Tx2>
	void projectionJL(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		constexpr Ti k_ = k * (1 + diffProj);

		Ti m = A.m, n = A.n;
		Tx T0 = Tx(0.0), T1 = Tx(1.0);
		Tx* d = o.d.get(), * L_d = o.L_d.get(), * AtL_d = o.AtL_d.get(), * x = o.x.get();

		for (Ti i = 0; i < m * k_; i++)
			d[i] = sign<Tx>(o.gen);

		for (Ti i = 0; i < n * k_; i++)
			AtL_d[i] = T0;

		ltsolve(L, (BaseImpl<Tx, k_>*)d, (BaseImpl<Tx, k_>*)L_d);
		gaxpy(At, (BaseImpl<Tx, k_>*)L_d, (BaseImpl<Tx, k_>*)AtL_d);

		for (Ti i = 0; i < n; i++)
		{
			Tx ret_i = T0;
			for (Ti j = 0; j < k; j++)
				ret_i += AtL_d[i * k_ + j] * AtL_d[i * k_ + j + k * diffProj];

			x[i] += ret_i;
		}
	}

	template<bool diffProj, typename Tx, typename Ti, typename Tx2>
	void projectionJL(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, size_t k)
	{
		if (!o.initialized())
			o.initialize(L, A, At);

		Ti n = A.n; Tx* x = o.x.get();
		for (Ti i = 0; i < n; i++)
			x[i] = Tx(0.0);

		constexpr size_t k_step = JLPackedSize / (1 + diffProj);
		for(size_t i = 1; i <= k / k_step; i++)
			projectionJL<diffProj, k_step>(o, L, A, At);

		for (size_t i = 1; i <= k % k_step; i++)
			projectionJL<diffProj, 1>(o, L, A, At);
	}

	template<typename Tx, typename Ti, typename Tx2>
	void leverageJL(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, size_t k)
	{
		projectionJL<false>(o, L, A, At, k);

		Ti n = o.n; Tx* x = o.x.get();
		Tx ratio = Tx(1 / double(k));
		for (Ti i = 0; i < n; i++)
			x[i] *= ratio;
	}

	template<typename Tx, typename Ti, typename Tx2>
	LeverageJLOutput<Tx, Ti> leverageJL(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, size_t k)
	{
		LeverageJLOutput<Tx, Ti> o;
		leverageJL(o, L, A, At, k);
		return o;
	}

	template<typename Tx, typename Ti, typename Tx2>
	Tx cholAccuracy(LeverageJLOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At, const Tx* w, size_t k = 1)
	{
		projectionJL<true>(o, L, A, At, k);

		Ti m = A.m, n = A.n;
		Tx* x = o.x.get();

		Tx ret1 = Tx(0.0);
		for (Ti i = 0; i < n; i++)
			ret1 += x[i] * w[i];

		Tx ret2 = Tx(0.0);
		for (Ti i = 0; i < m; i++)
		{
			Tx* d = o.d.get() + i * (2 * k);
			for (Ti j = 0; j < k; j++)
				ret2 += d[j] * d[j + k];

		}

		Tx dist = (ret1 - ret2) / Tx(sqrt(double(k)));
		return abs(dist);
	}
}