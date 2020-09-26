// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_USOLVE_H_
#define CSPARSE_CS_USOLVE_H_
#include "cs.h"

/* solve Ux=b where x and b are dense.  x=b on input, solution on output. */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_usolve(const cs<Tv,Ti>* U, Tv* x)
	{
		Ti p, j, n, * Up, * Ui;
		Tv* Ux;
		if (!CS_CSC(U) || !x) return (0);                     /* check inputs */
		n = U->n; Up = U->p; Ui = U->i; Ux = U->x;
		for (j = n - 1; j >= 0; j--)
		{
			x[j] /= Ux[Up[j + 1] - 1];
			for (p = Up[j]; p < Up[j + 1] - 1; p++)
			{
				x[Ui[p]] -= Ux[p] * x[j];
			}
		}
		return (1);
	}
}
#endif