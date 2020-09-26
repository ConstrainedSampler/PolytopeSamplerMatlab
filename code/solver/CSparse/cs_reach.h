// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_REACH_H_
#define CSPARSE_CS_REACH_H_
#include "cs.h"
#include "cs_dfs.h"

/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
namespace CSparse {
	template <typename Tv, typename Ti>
	Ti cs_reach(cs<Tv,Ti>* G, const cs<Tv, Ti>* B, Ti k, Ti* xi, const Ti* pinv)
	{
		Ti p, n, top, * Bp, * Bi, * Gp;
		if (!CS_CSC(G) || !CS_CSC(B) || !xi) return (-1);    /* check inputs */
		n = G->n; Bp = B->p; Bi = B->i; Gp = G->p;
		top = n;
		for (p = Bp[k]; p < Bp[k + 1]; p++)
		{
			if (!CS_MARKED(Gp, Bi[p]))    /* start a dfs at unmarked node i */
			{
				top = cs_dfs(Bi[p], G, top, xi, xi + n, pinv);
			}
		}
		for (p = top; p < n; p++) CS_MARK(Gp, xi[p]);  /* restore G */
		return (top);
	}
}
#endif
