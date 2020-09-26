// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SCATTER_H_
#define CSPARSE_CS_SCATTER_H_
#include "cs.h"

namespace CSparse {
	/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
	template <typename Tv, typename Ti>
	Ti cs_scatter(const cs<Tv,Ti>* A, Ti j, Tv beta, Ti* w, Tv* x, Ti mark,
		cs<Tv, Ti>* C, Ti nz)
	{
		Ti i, p, * Ap, * Ai, * Ci;
		Tv* Ax;
		if (!CS_CSC(A) || !w || !CS_CSC(C)) return (-1);     /* check inputs */
		Ap = A->p; Ai = A->i; Ax = A->x; Ci = C->i;
		for (p = Ap[j]; p < Ap[j + 1]; p++)
		{
			i = Ai[p];                            /* A(i,j) is nonzero */
			if (w[i] < mark)
			{
				w[i] = mark;                      /* i is new entry in column j */
				Ci[nz++] = i;                     /* add i to pattern of C(:,j) */
				if (x) x[i] = beta * Ax[p];      /* x(i) = beta*A(i,j) */
			}
			else if (x) x[i] += beta * Ax[p];    /* i exists in C(:,j) already */
		}
		return (nz);
	}
}

#endif