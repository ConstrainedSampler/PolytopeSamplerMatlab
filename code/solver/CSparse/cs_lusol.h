// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_LUSOL_H_
#define CSPARSE_CS_LUSOL_H_
#include "cs.h"
#include "cs_ipvec.h"
#include "cs_lsolve.h"
#include "cs_usolve.h"
#include "cs_sqr.h"
#include "cs_lu.h"

/* x=A\b where A is unsymmetric; b overwritten with solution */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_lusol(bool order, const cs<Tv,Ti>* A, Tv* b, Tv tol)
	{
		Tv* x;
		css<Tv, Ti>* S;
		csn<Tv, Ti>* N;
		Ti n;
		bool ok;
		if (!CS_CSC(A) || !b) return (0);     /* check inputs */
		n = A->n;
		S = cs_sqr(order, A, 0);              /* ordering and symbolic analysis */
		N = cs_lu(A, S, tol);                 /* numeric LU factorization */
		x = cs_malloc<Tv>(n);    /* get workspace */
		ok = (S && N && x);
		if (ok)
		{
			cs_ipvec(N->pinv, b, x, n);       /* x = b(p) */
			cs_lsolve(N->L, x);               /* x = L\x */
			cs_usolve(N->U, x);               /* x = U\x */
			cs_ipvec(S->q, x, b, n);          /* b(q) = x */
		}
		cs_free(x);
		cs_sfree(S);
		cs_nfree(N);
		return (ok);
	}
}
#endif