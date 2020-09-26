// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_SCHOL_H_
#define CSPARSE_CS_SCHOL_H_
#include "cs.h"
#include "cs_amd.h"
#include "cs_cumsum.h"
#include "cs_symperm.h"
#include "cs_pinv.h"
#include "cs_etree.h"
#include "cs_counts.h"
#include "cs_post.h"

/* ordering and symbolic analysis for a Cholesky factorization */
namespace CSparse {
	template <typename Tv, typename Ti>
	css<Tv,Ti>* cs_schol(bool order, const cs<Tv, Ti>* A)
	{
		Ti n, * c, * post, * P;
		cs<Tv, Ti>* C;
		css<Tv, Ti>* S;
		if (!CS_CSC(A)) return (NULL);        /* check inputs */
		n = A->n;
		S = cs_calloc<css<Tv, Ti>>(1);       /* allocate result S */
		if (!S) return (NULL);                 /* out of memory */
		P = cs_amd<Tv,Ti>(order, A);                 /* P = amd(A+A'), or natural */
		S->pinv = cs_pinv(P, n);              /* find inverse permutation */
		cs_free(P);
		if (order && !S->pinv) return (cs_sfree(S));
		C = cs_symperm(A, S->pinv, 0);        /* C = spones(triu(A(P,P))) */
		S->parent = cs_etree(C, 0);           /* find etree of C */
		post = cs_post(S->parent, n);         /* postorder the etree */
		c = cs_counts(C, S->parent, post, 0); /* find column counts of chol(C) */
		cs_free(post);
		cs_spfree(C);
		S->cp = cs_malloc<Ti>(n + 1); /* allocate result S->cp */
		S->unz = S->lnz = cs_cumsum<Ti>(S->cp, c, n); /* find column pointers for L */
		cs_free(c);
		return ((S->lnz >= 0) ? S : cs_sfree(S));
	}
}
#endif