// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SYMPERM_H_
#define CSPARSE_CS_SYMPERM_H_
#include "cs.h"
#include "cs_cumsum.h"

/* C = A(p,p) where A and C are symmetric the upper part stored; pinv not p */
namespace CSparse {
	template <typename Tv, typename Ti>
	cs<Tv,Ti>* cs_symperm(const cs<Tv, Ti>* A, const Ti* pinv, bool values)
	{
		Ti i, j, p, q, i2, j2, n, * Ap, * Ai, * Cp, * Ci, * w;
		Tv* Cx, * Ax;
		cs<Tv, Ti>* C;
		if (!CS_CSC(A)) return (NULL);                    /* check inputs */
		n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		C = cs_spalloc<Tv,Ti>(n, n, Ap[n], values && (Ax != NULL), 0); /* alloc result*/
		w = cs_calloc<Ti>(n);                   /* get workspace */
		if (!C || !w) return (cs_done(C, w, NULL, 0));    /* out of memory */
		Cp = C->p; Ci = C->i; Cx = C->x;
		for (j = 0; j < n; j++)           /* count entries in each column of C */
		{
			j2 = pinv ? pinv[j] : j;      /* column j of A is column j2 of C */
			for (p = Ap[j]; p < Ap[j + 1]; p++)
			{
				i = Ai[p];
				if (i > j) continue;       /* skip lower triangular part of A */
				i2 = pinv ? pinv[i] : i;  /* row i of A is row i2 of C */
				w[CS_MAX(i2, j2)]++;     /* column count of C */
			}
		}
		cs_cumsum<Ti>(Cp, w, n);              /* compute column pointers of C */
		for (j = 0; j < n; j++)
		{
			j2 = pinv ? pinv[j] : j;      /* column j of A is column j2 of C */
			for (p = Ap[j]; p < Ap[j + 1]; p++)
			{
				i = Ai[p];
				if (i > j) continue;       /* skip lower triangular part of A*/
				i2 = pinv ? pinv[i] : i;  /* row i of A is row i2 of C */
				Ci[q = w[CS_MAX(i2, j2)]++] = CS_MIN(i2, j2);
				if (Cx) Cx[q] = Ax[p];
			}
		}
		return (cs_done(C, w, NULL, 1));  /* success; free workspace, return C */
	}
}
#endif