// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_PERMUTE_H_
#define CSPARSE_CS_PERMUTE_H_
#include "cs.h"

/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1. */
namespace CSparse {
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_permute(const cs<Tv,Ti>* A, const Ti* pinv, const Ti* q, bool values)
	{
		Ti t, j, k, nz = 0, m, n, * Ap, * Ai, * Cp, * Ci;
		Tv* Cx, * Ax;
		cs<Tv, Ti>* C;
		if (!CS_CSC(A)) return (NULL);    /* check inputs */
		m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		C = cs_spalloc<Tv,Ti>(m, n, Ap[n], values && Ax != NULL, 0);  /* alloc result */
		if (!C) return (cs_done(C, NULL, NULL, 0));   /* out of memory */
		Cp = C->p; Ci = C->i; Cx = C->x;
		for (k = 0; k < n; k++)
		{
			Cp[k] = nz;                   /* column k of C is column q[k] of A */
			j = q ? (q[k]) : k;
			for (t = Ap[j]; t < Ap[j + 1]; t++)
			{
				if (Cx) Cx[nz] = Ax[t];  /* row i of A is row pinv[i] of C */
				Ci[nz++] = pinv ? (pinv[Ai[t]]) : Ai[t];
			}
		}
		Cp[n] = nz;                       /* finalize the last column of C */
		return (cs_done(C, NULL, NULL, 1));
	}
}
#endif
