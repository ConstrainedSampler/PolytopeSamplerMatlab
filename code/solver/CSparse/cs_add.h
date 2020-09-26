// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_ADD_H_
#define CSPARSE_CS_ADD_H_
#include "cs.h"
#include "cs_scatter.h"

namespace CSparse {
	/* C = alpha*A + beta*B */
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_add(const cs<Tv, Ti>* A, const cs<Tv, Ti>* B, Tv alpha, Tv beta)
	{
		Ti p, j, nz = 0, anz, * Cp, * Ci, * Bp, m, n, bnz, * w, values;
		Tv* x, * Bx, * Cx;
		cs<Tv, Ti>* C;
		if (!CS_CSC(A) || !CS_CSC(B)) return (NULL);         /* check inputs */
		if (A->m != B->m || A->n != B->n) return (NULL);
		m = A->m; anz = A->p[A->n];
		n = B->n; Bp = B->p; Bx = B->x; bnz = Bp[n];
		w = cs_calloc<Ti>(m);                       /* get workspace */
		values = (A->x != NULL) && (Bx != NULL);
		x = values ? cs_malloc<Tv>(m) : NULL;    /* get workspace */
		C = cs_spalloc<Tv, Ti>(m, n, anz + bnz, values, 0);           /* allocate result*/
		if (!C || !w || (values && !x)) return (cs_done<Tv, Ti>(C, w, x, 0));
		Cp = C->p; Ci = C->i; Cx = C->x;
		for (j = 0; j < n; j++)
		{
			Cp[j] = nz;                   /* column j of C starts here */
			nz = cs_scatter(A, j, alpha, w, x, j + 1, C, nz);   /* alpha*A(:,j)*/
			nz = cs_scatter(B, j, beta, w, x, j + 1, C, nz);    /* beta*B(:,j) */
			if (values) for (p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
		}
		Cp[n] = nz;                       /* finalize the last column of C */
		cs_sprealloc<Tv, Ti>(C, 0);               /* remove extra space from C */
		return (cs_done<Tv, Ti>(C, w, x, 1));     /* success; free workspace, return C */
	}
}

#endif