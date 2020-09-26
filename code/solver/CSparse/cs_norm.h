// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_NORM_H_
#define CSPARSE_CS_NORM_H_
#include "cs.h"

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
namespace CSparse {
	template <typename Tv, typename Ti>
	Tv cs_norm(const cs<Tv,Ti>* A)
	{
		Ti p, j, n, * Ap;
		Tv* Ax, norm = 0, s;
		if (!CS_CSC(A) || !A->x) return (-1);             /* check inputs */
		n = A->n; Ap = A->p; Ax = A->x;
		for (j = 0; j < n; j++)
		{
			for (s = 0, p = Ap[j]; p < Ap[j + 1]; p++) s += fabs(Ax[p]);
			norm = CS_MAX(norm, s);
		}
		return (norm);
	}
}
#endif