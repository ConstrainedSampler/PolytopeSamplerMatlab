// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_QRSOL_H_
#define CSPARSE_CS_QRSOL_H_
#include "cs.h"
#include "cs_happly.h"
#include "cs_sqr.h"
#include "cs_qr.h"
#include "cs_usolve.h"
#include "cs_utsolve.h"
#include "cs_ipvec.h"
#include "cs_pvec.h"
#include "cs_transpose.h"

/* x=A\b where A can be rectangular; b overwritten with solution */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_qrsol(bool order, const cs<Tv, Ti>* A, Tv* b)
	{
		Tv* x;
		css<Tv, Ti>* S;
		csn<Tv, Ti>* N;
		cs<Tv, Ti>* AT = NULL;
		Ti k, m, n, ok;
		if (!CS_CSC(A) || !b) return (0); /* check inputs */
		n = A->n;
		m = A->m;
		if (m >= n)
		{
			S = cs_sqr(order, A, 1);          /* ordering and symbolic analysis */
			N = cs_qr(A, S);                  /* numeric QR factorization */
			x = cs_calloc<Tv>(S ? S->m2 : 1);    /* get workspace */
			ok = (S && N && x);
			if (ok)
			{
				cs_ipvec(S->pinv, b, x, m);   /* x(0:m-1) = b(p(0:m-1) */
				for (k = 0; k < n; k++)       /* apply Householder refl. to x */
				{
					cs_happly(N->L, k, N->B[k], x);
				}
				cs_usolve(N->U, x);           /* x = R\x */
				cs_ipvec(S->q, x, b, n);      /* b(q(0:n-1)) = x(0:n-1) */
			}
		}
		else
		{
			AT = cs_transpose(A, 1);          /* Ax=b is underdetermined */
			S = cs_sqr(order, AT, 1);         /* ordering and symbolic analysis */
			N = cs_qr(AT, S);                 /* numeric QR factorization of A' */
			x = cs_calloc<Tv>(S ? S->m2 : 1);    /* get workspace */
			ok = (AT && S && N && x);
			if (ok)
			{
				cs_pvec(S->q, b, x, m);       /* x(q(0:m-1)) = b(0:m-1) */
				cs_utsolve(N->U, x);          /* x = R'\x */
				for (k = m - 1; k >= 0; k--)    /* apply Householder refl. to x */
				{
					cs_happly(N->L, k, N->B[k], x);
				}
				cs_pvec(S->pinv, x, b, n);    /* b(0:n-1) = x(p(0:n-1)) */
			}
		}
		cs_free(x);
		cs_sfree(S);
		cs_nfree(N);
		cs_spfree(AT);
		return (ok);
	}
}
#endif