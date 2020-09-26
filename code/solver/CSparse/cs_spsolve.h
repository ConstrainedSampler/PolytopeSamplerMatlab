// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SPSOLVE_H_
#define CSPARSE_CS_SPSOLVE_H_
#include "cs.h"
#include "cs_reach.h"

/* solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular */
namespace CSparse {
	template <typename Tv, typename Ti>
	Ti cs_spsolve(cs<Tv,Ti>* G, const cs<Tv, Ti>* B, Ti k, Ti* xi, Tv* x, const Ti* pinv, bool lo)
	{
		Ti j, J, p, q, px, top, n, * Gp, * Gi, * Bp, * Bi;
		Tv* Gx, * Bx;
		if (!CS_CSC(G) || !CS_CSC(B) || !xi || !x) return (-1);
		Gp = G->p; Gi = G->i; Gx = G->x; n = G->n;
		Bp = B->p; Bi = B->i; Bx = B->x;
		top = cs_reach(G, B, k, xi, pinv);        /* xi[top..n-1]=Reach(B(:,k)) */
		for (p = top; p < n; p++) x[xi[p]] = Tv(0);    /* clear x */
		for (p = Bp[k]; p < Bp[k + 1]; p++) x[Bi[p]] = Bx[p]; /* scatter B */
		for (px = top; px < n; px++)
		{
			j = xi[px];                               /* x(j) is nonzero */
			J = pinv ? (pinv[j]) : j;                 /* j maps to col J of G */
			if (J < 0) continue;                       /* column J is empty */
			x[j] /= Gx[lo ? (Gp[J]) : (Gp[J + 1] - 1)];/* x(j) /= G(j,j) */
			p = lo ? (Gp[J] + 1) : (Gp[J]);            /* lo: L(j,j) 1st entry */
			q = lo ? (Gp[J + 1]) : (Gp[J + 1] - 1);        /* up: U(j,j) last entry */
			for (; p < q; p++)
			{
				x[Gi[p]] -= Gx[p] * x[j];          /* x(i) -= G(i,j) * x(j) */
			}
		}
		return (top);                                  /* return top of stack */
	}
}
#endif
