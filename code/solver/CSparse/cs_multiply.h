// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_MULTIPLY_H_
#define CSPARSE_CS_MULTIPLY_H_
#include "cs.h"
#include "cs_scatter.h"

/* C = A*B */
namespace CSparse {
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_multiply(const cs<Tv, Ti>* A, const cs<Tv, Ti>* B, bool values = true)
	{
		Ti p, j, nz = 0, anz, * Cp, * Ci, * Bp, m, n, bnz, * w, * Bi;
		Tv* x, * Bx, * Cx;
		cs<Tv, Ti>* C;
		if (!CS_CSC(A) || !CS_CSC(B)) return (NULL);      /* check inputs */
		if (A->n != B->m) return (NULL);
		m = A->m; anz = A->p[A->n];
		n = B->n; Bp = B->p; Bi = B->i; Bx = B->x; bnz = Bp[n];
		w = cs_calloc<Ti>(m);                    /* get workspace */
		values = (A->x != NULL) && (Bx != NULL) && values;
		x = values ? cs_malloc<Tv>(m) : NULL; /* get workspace */
		C = cs_spalloc<Tv,Ti>(m, n, anz + bnz, values, 0);        /* allocate result */
		if (!C || !w || (values && !x)) return (cs_done(C, w, x, 0));
		Cp = C->p;
		for (j = 0; j < n; j++)
		{
			if (nz + m > C->nzmax && !cs_sprealloc(C, 2 * (C->nzmax) + m))
			{
				return (cs_done(C, w, x, 0));             /* out of memory */
			}
			Ci = C->i; Cx = C->x;         /* C->i and C->x may be reallocated */
			Cp[j] = nz;                   /* column j of C starts here */
			for (p = Bp[j]; p < Bp[j + 1]; p++)
			{
				nz = cs_scatter(A, Bi[p], Bx ? Bx[p] : 1, w, x, j + 1, C, nz);
			}
			if (values) for (p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
		}
		Cp[n] = nz;                       /* finalize the last column of C */
		cs_sprealloc<Tv,Ti>(C, 0);               /* remove extra space from C */
		return (cs_done(C, w, x, 1));     /* success; free workspace, return C */
	}
}
#endif