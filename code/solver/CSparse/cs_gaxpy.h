// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_GAXPY_H_
#define CSPARSE_CS_GAXPY_H_
#include "cs.h"
/* y = A*x+y */

namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_gaxpy(const cs<Tv,Ti>* A, const Tv* x, Tv* y)
	{
		Ti p, j, n, * Ap, * Ai;
		Tv* Ax;
		if (!CS_CSC(A) || !x || !y) return (0);       /* check inputs */
		n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		for (j = 0; j < n; j++)
		{
			for (p = Ap[j]; p < Ap[j + 1]; p++)
			{
				y[Ai[p]] += Ax[p] * x[j];
			}
		}
		return (1);
	}
}
#endif
