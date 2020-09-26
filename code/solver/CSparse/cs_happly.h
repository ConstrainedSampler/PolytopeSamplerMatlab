// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_HAPPLY_H_
#define CSPARSE_CS_HAPPLY_H_
#include "cs.h"

/* apply the ith Householder vector to x */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_happly(const cs<Tv,Ti>* V, Ti i, Tv beta, Tv* x)
	{
		Ti p, * Vp, * Vi;
		Tv* Vx, tau = 0;
		if (!CS_CSC(V) || !x) return (0);     /* check inputs */
		Vp = V->p; Vi = V->i; Vx = V->x;
		for (p = Vp[i]; p < Vp[i + 1]; p++)   /* tau = v'*x */
		{
			tau += Vx[p] * x[Vi[p]];
		}
		tau *= beta;                           /* tau = beta*(v'*x) */
		for (p = Vp[i]; p < Vp[i + 1]; p++)   /* x = x - v*tau */
		{
			x[Vi[p]] -= Vx[p] * tau;
		}
		return (1);
	}
}
#endif