// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_PRINT_H_
#define CSPARSE_CS_PRINT_H_
#include "cs.h"
#include "cs_norm.h"
#include <cstdio>

/* print a sparse matrix; use %g for integers to avoid differences with csi */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_print(const cs<Tv,Ti>* A, bool brief)
	{
		Ti p, j, m, n, nzmax, * Ap, * Ai;
		Tv* Ax;
		if (!A) { printf("(null)\n"); return (0); }
		m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		nzmax = A->nzmax;
		printf("CSparse Version %d.%d.%d, %s.  %s\n", 3, 2, 0, "Sept 12, 2017", "Copyright (c) Timothy A. Davis, 2006-2016");
		printf("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double)m,
			(double)n, (double)nzmax, (double)(Ap[n]), cs_norm(A));
		for (j = 0; j < n; j++)
		{
			printf("    col %g : locations %g to %g\n", (double)j,
				(double)(Ap[j]), (double)(Ap[j + 1] - 1));
			for (p = Ap[j]; p < Ap[j + 1]; p++)
			{
				printf("      %g : %g\n", (double)(Ai[p]), Ax ? Ax[p] : 1);
				if (brief && p > 20) { printf("  ...\n"); return (1); }
			}
		}
		return (1);
	}
}

#endif