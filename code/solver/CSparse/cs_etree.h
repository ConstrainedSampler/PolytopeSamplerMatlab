// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_ETREE_H_
#define CSPARSE_CS_ETREE_H_
#include "cs.h"

/* compute the etree of A (using triu(A), or A'A without forming A'A */
namespace CSparse {
	template <typename Tv, typename Ti>
	Ti* cs_etree(const cs<Tv,Ti>* A, bool ata)
	{
		Ti i, k, p, m, n, inext, * Ap, * Ai, * w, * parent, * ancestor, * prev;
		if (!CS_CSC(A)) return (NULL);        /* check inputs */
		m = A->m; n = A->n; Ap = A->p; Ai = A->i;
		parent = cs_malloc<Ti>(n);              /* allocate result */
		w = cs_malloc<Ti>(n + (ata ? m : 0));   /* get workspace */
		if (!w || !parent) return (cs_idone<Tv,Ti>(parent, NULL, w, 0));
		ancestor = w; prev = w + n;
		if (ata) for (i = 0; i < m; i++) prev[i] = -1;
		for (k = 0; k < n; k++)
		{
			parent[k] = -1;                   /* node k has no parent yet */
			ancestor[k] = -1;                 /* nor does k have an ancestor */
			for (p = Ap[k]; p < Ap[k + 1]; p++)
			{
				i = ata ? (prev[Ai[p]]) : (Ai[p]);
				for (; i != -1 && i < k; i = inext)   /* traverse from i to k */
				{
					inext = ancestor[i];              /* inext = ancestor of i */
					ancestor[i] = k;                  /* path compression */
					if (inext == -1) parent[i] = k;   /* no anc., parent is k */
				}
				if (ata) prev[Ai[p]] = k;
			}
		}
		return (cs_idone<Tv,Ti>(parent, NULL, w, 1));
	}
}

#endif