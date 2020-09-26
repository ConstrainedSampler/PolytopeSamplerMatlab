// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SCC_H_
#define CSPARSE_CS_SCC_H_
#include "cs.h"
#include "cs_dfs.h"

/* find the strongly connected components of a square matrix */
namespace CSparse {
	template <typename Tv, typename Ti>
	csd<Ti>* cs_scc(cs<Tv,Ti>* A)     /* matrix A temporarily modified, then restored */
	{
		Ti n, i, k, b, nb = 0, top, * xi, * pstack, * p, * r, * Ap, * ATp, * rcopy, * Blk;
		cs<Tv, Ti>* AT;
		csd<Ti>* D;
		if (!CS_CSC(A)) return (NULL);                /* check inputs */
		n = A->n; Ap = A->p;
		D = cs_dalloc<Ti>(n, 0);                          /* allocate result */
		AT = cs_transpose(A, 0);                      /* AT = A' */
		xi = cs_malloc<Ti>(2 * n + 1);          /* get workspace */
		if (!D || !AT || !xi) return (cs_ddone(D, AT, xi, 0));
		Blk = xi; rcopy = pstack = xi + n;
		p = D->p; r = D->r; ATp = AT->p;
		top = n;
		for (i = 0; i < n; i++)   /* first dfs(A) to find finish times (xi) */
		{
			if (!CS_MARKED(Ap, i)) top = cs_dfs<Tv,Ti>(i, A, top, xi, pstack, NULL);
		}
		for (i = 0; i < n; i++) CS_MARK(Ap, i); /* restore A; unmark all nodes*/
		top = n;
		nb = n;
		for (k = 0; k < n; k++)   /* dfs(A') to find strongly connnected comp */
		{
			i = xi[k];            /* get i in reverse order of finish times */
			if (CS_MARKED(ATp, i)) continue;  /* skip node i if already ordered */
			r[nb--] = top;        /* node i is the start of a component in p */
			top = cs_dfs<Tv,Ti>(i, AT, top, p, pstack, NULL);
		}
		r[nb] = 0;                /* first block starts at zero; shift r up */
		for (k = nb; k <= n; k++) r[k - nb] = r[k];
		D->nb = nb = n - nb;         /* nb = # of strongly connected components */
		for (b = 0; b < nb; b++)  /* sort each block in natural order */
		{
			for (k = r[b]; k < r[b + 1]; k++) Blk[p[k]] = b;
		}
		for (b = 0; b <= nb; b++) rcopy[b] = r[b];
		for (i = 0; i < n; i++) p[rcopy[Blk[i]]++] = i;
		return (cs_ddone(D, AT, xi, 1));
	}
}
#endif