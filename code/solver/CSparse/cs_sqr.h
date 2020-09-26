// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SQR_H_
#define CSPARSE_CS_SQR_H_
#include "cs.h"
#include "cs_amd.h"
#include "cs_etree.h"
#include "cs_permute.h"
#include "cs_post.h"
#include "cs_counts.h"


namespace CSparse {
	/* compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent */
	template <typename Tv, typename Ti>
	bool cs_vcount(const cs<Tv,Ti>* A, css<Tv, Ti>* S)
	{
		Ti i, k, p, pa, n = A->n, m = A->m, * Ap = A->p, * Ai = A->i, * next, * head,
			* tail, * nque, * pinv, * leftmost, * w, * parent = S->parent;
		S->pinv = pinv = cs_malloc<Ti>(m + n);        /* allocate pinv, */
		S->leftmost = leftmost = cs_malloc<Ti>(m);  /* and leftmost */
		w = cs_malloc<Ti>(m + 3 * n);   /* get workspace */
		if (!pinv || !w || !leftmost)
		{
			cs_free(w);                       /* pinv and leftmost freed later */
			return (0);                        /* out of memory */
		}
		next = w; head = w + m; tail = w + m + n; nque = w + m + 2 * n;
		for (k = 0; k < n; k++) head[k] = -1;   /* queue k is empty */
		for (k = 0; k < n; k++) tail[k] = -1;
		for (k = 0; k < n; k++) nque[k] = 0;
		for (i = 0; i < m; i++) leftmost[i] = -1;
		for (k = n - 1; k >= 0; k--)
		{
			for (p = Ap[k]; p < Ap[k + 1]; p++)
			{
				leftmost[Ai[p]] = k;         /* leftmost[i] = min(find(A(i,:)))*/
			}
		}
		for (i = m - 1; i >= 0; i--)            /* scan rows in reverse order */
		{
			pinv[i] = -1;                     /* row i is not yet ordered */
			k = leftmost[i];
			if (k == -1) continue;             /* row i is empty */
			if (nque[k]++ == 0) tail[k] = i; /* first row in queue k */
			next[i] = head[k];               /* put i at head of queue k */
			head[k] = i;
		}
		S->lnz = 0;
		S->m2 = m;
		for (k = 0; k < n; k++)               /* find row permutation and nnz(V)*/
		{
			i = head[k];                      /* remove row i from queue k */
			S->lnz++;                          /* count V(k,k) as nonzero */
			if (i < 0) i = S->m2++;            /* add a fictitious row */
			pinv[i] = k;                      /* associate row i with V(:,k) */
			if (--nque[k] <= 0) continue;     /* skip if V(k+1:m,k) is empty */
			S->lnz += nque[k];                /* nque [k] is nnz (V(k+1:m,k)) */
			if ((pa = parent[k]) != -1)        /* move all rows to parent of k */
			{
				if (nque[pa] == 0) tail[pa] = tail[k];
				next[tail[k]] = head[pa];
				head[pa] = next[i];
				nque[pa] += nque[k];
			}
		}
		for (i = 0; i < m; i++) if (pinv[i] < 0) pinv[i] = k++;
		cs_free(w);
		return (1);
	}

	/* symbolic ordering and analysis for QR or LU */
	template <typename Tv, typename Ti>
	css<Tv, Ti>* cs_sqr(bool order, const cs<Tv, Ti>* A, bool qr)
	{
		Ti n, k, ok = 1, * post;
		css<Tv, Ti>* S;
		if (!CS_CSC(A)) return (NULL);        /* check inputs */
		n = A->n;
		S = cs_calloc<css<Tv, Ti>>(1);       /* allocate result S */
		if (!S) return (NULL);                 /* out of memory */
		S->q = cs_amd<Tv,Ti>(order, A);              /* fill-reducing ordering */
		if (order && !S->q) return (cs_sfree(S));
		if (qr)                                 /* QR symbolic analysis */
		{
			cs<Tv, Ti>* C = order ? cs_permute<Tv,Ti>(A, NULL, S->q, 0) : ((cs<Tv,Ti>*)A);
			S->parent = cs_etree<Tv,Ti>(C, 1);       /* etree of C'*C, where C=A(:,q) */
			post = cs_post(S->parent, n);
			S->cp = cs_counts(C, S->parent, post, 1);  /* col counts chol(C'*C) */
			cs_free(post);
			ok = C && S->parent && S->cp && cs_vcount(C, S);
			if (ok) for (S->unz = 0, k = 0; k < n; k++) S->unz += S->cp[k];
			if (order) cs_spfree<Tv,Ti>(C);
		}
		else
		{
			S->unz = (int64_t)(4 * (A->p[n]) + n);         /* for LU factorization only, */
			S->lnz = (int64_t)S->unz;                   /* guess nnz(L) and nnz(U) */
		}
		return (ok ? S : cs_sfree(S));        /* return result S */
	}
}

#endif