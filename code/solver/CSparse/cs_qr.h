// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_QR_H_
#define CSPARSE_CS_QR_H_
#include "cs.h"
#include "cs_scatter.h"
#include "cs_house.h"

/* sparse QR factorization [V,beta,pinv,R] = qr (A) */
namespace CSparse {
	template <typename Tv, typename Ti>
	csn<Tv,Ti>* cs_qr(const cs<Tv, Ti>* A, const css<Tv, Ti>* S)
	{
		Tv* Rx, * Vx, * Ax, * x, * Beta;
		Ti i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, * s, * leftmost, * Ap, * Ai,
			* parent, * Rp, * Ri, * Vp, * Vi, * w, * pinv, * q;
		cs<Tv, Ti>* R, * V;
		csn<Tv, Ti>* N;
		if (!CS_CSC(A) || !S) return (NULL);
		m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		q = S->q; parent = S->parent; pinv = S->pinv; m2 = S->m2;
		vnz = (Ti)S->lnz; rnz = (Ti)S->unz; leftmost = S->leftmost;
		w = cs_malloc<Ti>(m2 + n);            /* get csi workspace */
		x = cs_malloc<Tv>(m2);           /* get double workspace */
		N = cs_calloc<csn<Tv,Ti>>(1);               /* allocate result */
		if (!w || !x || !N) return (cs_ndone<Tv,Ti>(N, NULL, w, x, 0));
		s = w + m2;                                    /* s is size n */
		for (k = 0; k < m2; k++) x[k] = Tv(0);          /* clear workspace x */
		N->L = V = cs_spalloc<Tv, Ti>(m2, n, vnz, 1, 0);      /* allocate result V */
		N->U = R = cs_spalloc<Tv, Ti>(m2, n, rnz, 1, 0);      /* allocate result R */
		N->B = Beta = cs_malloc<Tv>(n);  /* allocate result Beta */
		if (!R || !V || !Beta) return (cs_ndone<Tv, Ti>(N, NULL, w, x, 0));
		Rp = R->p; Ri = R->i; Rx = R->x;
		Vp = V->p; Vi = V->i; Vx = V->x;
		for (i = 0; i < m2; i++) w[i] = -1; /* clear w, to mark nodes */
		rnz = 0; vnz = 0;
		for (k = 0; k < n; k++)               /* compute V and R */
		{
			Rp[k] = rnz;                      /* R(:,k) starts here */
			Vp[k] = p1 = vnz;                 /* V(:,k) starts here */
			w[k] = k;                         /* add V(k,k) to pattern of V */
			Vi[vnz++] = k;
			top = n;
			col = q ? q[k] : k;
			for (p = Ap[col]; p < Ap[col + 1]; p++)   /* find R(:,k) pattern */
			{
				i = leftmost[Ai[p]];         /* i = min(find(A(i,q))) */
				for (len = 0; w[i] != k; i = parent[i]) /* traverse up to k */
				{
					s[len++] = i;
					w[i] = k;
				}
				while (len > 0) s[--top] = s[--len]; /* push path on stack */
				i = pinv[Ai[p]];             /* i = permuted row of A(:,col) */
				x[i] = Ax[p];                /* x (i) = A(:,col) */
				if (i > k&& w[i] < k)         /* pattern of V(:,k) = x (k+1:m) */
				{
					Vi[vnz++] = i;            /* add i to pattern of V(:,k) */
					w[i] = k;
				}
			}
			for (p = top; p < n; p++) /* for each i in pattern of R(:,k) */
			{
				i = s[p];                     /* R(i,k) is nonzero */
				cs_happly(V, i, Beta[i], x); /* apply (V(i),Beta(i)) to x */
				Ri[rnz] = i;                  /* R(i,k) = x(i) */
				Rx[rnz++] = x[i];
				x[i] = Tv(0);
				if (parent[i] == k) vnz = cs_scatter<Tv,Ti>(V, i, 0, w, NULL, k, V, vnz);
			}
			for (p = p1; p < vnz; p++)        /* gather V(:,k) = x */
			{
				Vx[p] = x[Vi[p]];
				x[Vi[p]] = Tv(0);
			}
			Ri[rnz] = k;                     /* R(k,k) = norm (x) */
			Rx[rnz++] = cs_house(Vx + p1, Beta + k, vnz - p1); /* [v,beta]=house(x) */
		}
		Rp[n] = rnz;                          /* finalize R */
		Vp[n] = vnz;                          /* finalize V */
		return (cs_ndone<Tv, Ti>(N, NULL, w, x, 1));  /* success */
	}
}
#endif
