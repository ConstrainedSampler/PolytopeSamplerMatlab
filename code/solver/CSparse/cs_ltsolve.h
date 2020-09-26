// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_LTSOLVE_H_
#define CSPARSE_CS_LTSOLVE_H_
#include "cs.h"

/* solve L'x=b where x and b are dense.  x=b on input, solution on output. */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_ltsolve(const cs<Tv,Ti>* L, Tv* x)
	{
		Ti p, j, n, * Lp, * Li;
		Tv* Lx;
		if (!CS_CSC(L) || !x) return (0);                     /* check inputs */
		n = L->n; Lp = L->p; Li = L->i; Lx = L->x;
		for (j = n - 1; j >= 0; j--)
		{
			for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
			{
				x[j] -= Lx[p] * x[Li[p]];
			}
			x[j] /= Lx[Lp[j]];
		}
		return (1);
	}
}
#endif