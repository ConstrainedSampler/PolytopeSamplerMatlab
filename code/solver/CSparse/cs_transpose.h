// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_TRANSPOSE_H_
#define CSPARSE_CS_TRANSPOSE_H_
#include "cs.h"
#include "cs_cumsum.h"

namespace CSparse {
	/* C = A' */
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_transpose(const cs<Tv, Ti>* A, bool values)
	{
		Ti p, q, j, * Cp, * Ci, n, m, * Ap, * Ai, * w;
		Tv* Cx, * Ax;
		cs<Tv, Ti>* C;
		if (!CS_CSC(A)) return (NULL);    /* check inputs */
		m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		C = cs_spalloc<Tv, Ti>(n, m, Ap[n], values && Ax, 0);       /* allocate result */
		w = cs_calloc<Ti>(m);                      /* get workspace */
		if (!C || !w) return (cs_done<Tv, Ti>(C, w, NULL, 0));       /* out of memory */
		Cp = C->p; Ci = C->i; Cx = C->x;
		for (p = 0; p < Ap[n]; p++) w[Ai[p]]++;          /* row counts */
		cs_cumsum<Ti>(Cp, w, m);                                 /* row pointers */
		for (j = 0; j < n; j++)
		{
			for (p = Ap[j]; p < Ap[j + 1]; p++)
			{
				Ci[q = w[Ai[p]]++] = j; /* place A(i,j) as entry C(j,i) */
				if (Cx) Cx[q] = Ax[p];
			}
		}
		return (cs_done<Tv, Ti>(C, w, NULL, 1));  /* success; free w and return C */
	}
}

#endif