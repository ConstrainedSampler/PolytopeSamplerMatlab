// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_H_
#define CSPARSE_CS_H_
#include <cstdint>
#include <cmath>

namespace CSparse {
	/* --- primary CSparse data structures ------------------------- */
	template <typename Tv, typename Ti>
	struct cs			/* matrix in compressed-column or triplet form */
	{
		Ti nzmax;		/* maximum number of entries */
		Ti m;			/* number of rows */
		Ti n;			/* number of columns */
		Ti *p;			/* column pointers (size n+1) or col indices (size nzmax) */
		Ti *i;			/* row indices, size nzmax */
		Tv *x;			/* numerical values, size nzmax */
	};

	///* --- secondary CSparse data structures ----------------------- */
	template <typename Tv, typename Ti>
	struct css			/* symbolic Cholesky, LU, or QR analysis */
	{
		Ti *pinv;		/* inverse row perm. for QR, fill red. perm for Chol */
		Ti *q;			/* fill-reducing column permutation for LU and QR */
		Ti *parent;		/* elimination tree for Cholesky and QR */
		Ti *cp;			/* column pointers for Cholesky, row counts for QR */
		Ti *leftmost;	/* leftmost[i] = min(find(A(i,:))), for QR */
		Ti m2;			/* # of rows for QR, after adding fictitious rows */
		int64_t lnz;			/* # entries in L for LU or Cholesky; in V for QR */
		int64_t unz;			/* # entries in U for LU; in R for QR */
	};

	template <typename Tv, typename Ti>
	struct csn			/* numeric Cholesky, LU, or QR factorization */
	{
	    cs<Tv,Ti> *L;	/* L for LU and Cholesky, V for QR */
	    cs<Tv,Ti> *U;	/* U for LU, R for QR, not used for Cholesky */
	    Ti *pinv;		/* partial pivoting for LU */
	    Tv *B;			/* beta [0..n-1] for QR */
	};

	template <typename Ti>
	struct csd			/* cs_dmperm or cs_scc output */
	{
		Ti *p;			/* size m, row permutation */
		Ti *q;			/* size n, column permutation */
		Ti *r;			/* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
		Ti *s;			/* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
		Ti nb;			/* # of blocks in fine dmperm decomposition */
		Ti rr[5];		/* coarse row decomposition */
		Ti cc[5];		/* coarse column decomposition */
	};

	#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
	#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
	#define CS_FLIP(i) (-(i)-2)
	#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
	#define CS_MARKED(w,j) (w [j] < 0)
	#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
	#define CS_CSC(A) (A)

	template <typename T, typename Ti>
	T* cs_malloc(Ti n) // no initialization
	{
		return new T[n];
	}
    
	template <typename T, typename Ti>
	T* cs_calloc(Ti n)
	{
		return new T[n]();
	}
    
	template <typename T>
	T* cs_free(T* p)
	{
		if (p) delete[] p;       /* free p if it is not already NULL */
		return (nullptr);         /* return NULL to simplify the use of cs_free */
	}
    
	template <typename T>
	T* cs_realloc(T* p, size_t n_org, size_t n, bool* ok)
	{
        T* pnew = cs_calloc<T>(n);
        if (pnew)
        {
            size_t k = CS_MIN(n, n_org);
            for (size_t i = 0; i < k; ++i)
                pnew[i] = p[i];
            cs_free(p);
        }
		*ok = (pnew != nullptr);                  /* realloc fails if pnew is NULL */
		return ((*ok) ? pnew : p);             /* return original p if failure */
	}
    
    template <typename Tv_, typename Ti_, typename Tv, typename Ti>
    cs<Tv, Ti>* cs_copy(const cs<Tv_, Ti_>* A_)
    {

        cs<Tv, Ti>* A = cs_spalloc<Tv, Ti>(A_->m, A_->n, A_->nzmax, A_->x != nullptr, false);

        for (size_t s = 0; s <= A_->n; ++s)
            A->p[s] = convert<Ti_, Ti>(A_->p[s]);

        for (size_t s = 0; s < A_->nzmax; ++s)
            A->i[s] = convert<Ti_, Ti>(A_->i[s]);

        if (A_->x)
        {
            for (size_t s = 0; s < A_->nzmax; ++s)
                A->x[s] = convert<Tv_, Tv>(A_->x[s]);
        }

        return A;
    }
    
	///* allocate a sparse matrix (triplet form or compressed-column form) */
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_spalloc(Ti m, Ti n, Ti nzmax, bool values, bool triplet)
	{
		cs<Tv, Ti>* A = cs_calloc<cs<Tv, Ti>>(1);    /* allocate the cs struct */
		if (!A) return (NULL);                 /* out of memory */
		A->m = m;                              /* define dimensions and nzmax */
		A->n = n;
		A->nzmax = nzmax = CS_MAX(nzmax, 1);
		A->p = cs_malloc<Ti>(triplet ? nzmax : n + 1);
		A->i = cs_malloc<Ti>(nzmax);
		A->x = values ? cs_malloc<Tv>(nzmax) : NULL;
		return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree(A) : A);
	}

	/* change the max # of entries sparse matrix */
	template <typename Tv, typename Ti>
	bool cs_sprealloc(cs<Tv, Ti>* A, Ti nzmax)
	{
		bool ok, oki, okj = 1, okx = 1;
		if (!A) return (0);
		if (nzmax <= 0) nzmax = A->p[A->n];
		nzmax = CS_MAX(nzmax, 1);
		A->i = cs_realloc<Ti>(A->i, A->nzmax, nzmax, &oki);
		if (A->x) A->x = cs_realloc<Tv>(A->x, A->nzmax, nzmax, &okx);
		ok = (oki && okj && okx);
		if (ok) A->nzmax = nzmax;
		return (ok);
	}

	///* free a sparse matrix */
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_spfree(cs<Tv, Ti>* A)
	{
		if (!A) return (NULL);     /* do nothing if A already NULL */
		cs_free(A->p);
		cs_free(A->i);
		cs_free(A->x);
		return cs_free(A);   /* free the cs struct and return NULL */
	}

	///* free a numeric factorization */
	template <typename Tv, typename Ti>
	csn<Tv,Ti>* cs_nfree(csn<Tv, Ti>* N)
	{
		if (!N) return (NULL);     /* do nothing if N already NULL */
		cs_spfree(N->L);
		cs_spfree(N->U);
		cs_free(N->pinv);
		cs_free(N->B);
		return cs_free(N);  /* free the csn struct and return NULL */
	}

	///* free a symbolic factorization */
	template <typename Tv, typename Ti>
	css<Tv, Ti>* cs_sfree(css<Tv,Ti>* S)
	{
		if (!S) return (NULL);     /* do nothing if S already NULL */
		cs_free(S->pinv);
		cs_free(S->q);
		cs_free(S->parent);
		cs_free(S->cp);
		cs_free(S->leftmost);
		return cs_free(S);  /* free the css struct and return NULL */
	}

	///* allocate a cs_dmperm or cs_scc result */
	template <typename Ti>
	csd<Ti>* cs_dalloc(Ti m, Ti n)
	{
		csd<Ti>* D;
		D = cs_calloc<csd<Ti>>(1);
		if (!D) return (NULL);
		D->p = cs_malloc<Ti>(m);
		D->r = cs_malloc<Ti>(m + 6);
		D->q = cs_malloc<Ti>(n);
		D->s = cs_malloc<Ti>(n + 6);
		return ((!D->p || !D->r || !D->q || !D->s) ? cs_dfree(D) : D);
	}

	///* free a cs_dmperm or cs_scc result */
	template <typename Ti>
	csd<Ti>* cs_dfree(csd<Ti>* D)
	{
		if (!D) return (NULL);     /* do nothing if D already NULL */
		cs_free(D->p);
		cs_free(D->q);
		cs_free(D->r);
		cs_free(D->s);
		return cs_free(D);  /* free the csd struct and return NULL */
	}

	/* free workspace and return a sparse matrix result */
	template <typename Tv, typename Ti>
	cs<Tv, Ti>* cs_done(cs<Tv, Ti>* C, void* w, void* x, bool ok)
	{
		cs_free(w);                       /* free workspace */
		cs_free(x);
		return (ok ? C : cs_spfree(C));   /* return result if OK, else free it */
	}

	///* free workspace and return csi array result */
	template <typename Tv, typename Ti>
	Ti* cs_idone(Ti* p, cs<Tv,Ti>* C, void* w, bool ok)
	{
		cs_spfree(C);                     /* free temporary matrix */
		cs_free(w);                       /* free workspace */
		return (ok ? p : (Ti*)cs_free(p)); /* return result, or free it */
	}

	///* free workspace and return a numeric factorization (Cholesky, LU, or QR) */
	template <typename Tv, typename Ti>
	csn<Tv, Ti>* cs_ndone(csn<Tv, Ti>* N, cs<Tv, Ti>* C, void* w, void* x, bool ok)
	{
		cs_spfree(C);                     /* free temporary matrix */
		cs_free(w);                       /* free workspace */
		cs_free(x);
		return (ok ? N : cs_nfree(N));    /* return result if OK, else free it */
	}

	///* free workspace and return a csd result */
	template <typename Tv, typename Ti>
	csd<Ti>* cs_ddone(csd<Ti>* D, cs<Tv, Ti>* C, void* w, bool ok)
	{
		cs_spfree(C);                     /* free temporary matrix */
		cs_free(w);                       /* free workspace */
		return (ok ? D : cs_dfree(D));    /* return result if OK, else free it */
	}
}
#endif
