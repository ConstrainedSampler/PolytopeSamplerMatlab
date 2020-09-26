// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_CHOLSOL_H_
#define CSPARSE_CS_CHOLSOL_H_
#include "cs.h"
#include "cs_schol.h"
#include "cs_chol.h"
#include "cs_lsolve.h"
#include "cs_ltsolve.h"
#include "cs_ipvec.h"
#include "cs_pvec.h"

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_cholsol(bool order, const cs<Tv,Ti>* A, Tv* b)
	{
		Tv* x;
		css<Tv, Ti>* S;
		csn<Tv, Ti>* N;
		Ti n;
		bool ok;
		if (!CS_CSC(A) || !b) return (0);     /* check inputs */
		n = A->n;
		S = cs_schol(order, A);               /* ordering and symbolic analysis */
		N = cs_chol(A, S);                    /* numeric Cholesky factorization */
		x = cs_malloc<Tv>(n);    /* get workspace */
		ok = (S && N && x);
		if (ok)
		{
			cs_ipvec(S->pinv, b, x, n);   /* x = P*b */
			cs_lsolve(N->L, x);           /* x = L\x */
			cs_ltsolve(N->L, x);          /* x = L'\x */
			cs_pvec(S->pinv, x, b, n);    /* b = P'*x */
		}
		cs_free(x);
		cs_sfree(S);
		cs_nfree(N);
		return (ok);
	}
}
#endif