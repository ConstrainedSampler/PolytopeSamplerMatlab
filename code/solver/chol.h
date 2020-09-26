#include "CSparse/cs_chol.h"
#include "CSparse/cs_schol.h"
#include "CSparse/cs_ereach.h"

// This program assume the diag of H is non-zero.
// It will output incorrect answer if it is not the case.
template<typename Tv, typename Ti>
struct CholObj
{
    cs<bool, Ti> *Lt = nullptr;
    Ti *A_diag_ = nullptr; // Ax[A_diag_[k]] is A_kk
    
    // workspaces
    Tv *x_ = nullptr; 
    Ti *c_ = nullptr; // c[i] = index the last nonzero on column i in the current L
    
    cs<Tv, Ti>* prepare(const cs<Tv, Ti>* A)
    {
        css<Tv, Ti>* S = cs_schol(false, A);
        Ti *cp = S->cp, *parent = S->parent; Ti n = A->n;
        
        c_ = cs_malloc<Ti>(2 * n);
        x_ = cs_malloc<Tv>(n);
        cs<Tv, Ti>* out = cs_spalloc<Tv, Ti>(n, n, cp[n], true, false);
        
        Tv *x = x_; // Loop over variables in structure is slower somehow
        Ti *c = c_, *s = c_ + n;

        Ti *Lp = out->p, *Li = out->i;
        for (Ti k = 0; k < n; k++)
            Lp[k] = c[k] = cp[k];

        for (Ti k = 0; k < n; k++)
        {
            Ti top = cs_ereach<Tv,Ti>(A, k, parent, s, c);
            for (; top < n; top++)
            {
                Li[c[s[top]]++] = k;
            }
            Li[c[k]++] = k;
        }
        Lp[n] = cp[n];
        cs_sfree(S);
        
        
        Lt = cs_transpose<bool, Ti>((cs<bool, Ti>*)out, false);
        
        Ti *Ap = A->p, *Ai = A->i;
        Ti* A_diag = cs_malloc<Ti>(n);
        for (Ti k = 0; k < n; ++k)
        {
            Ti s = Ap[k], s_end = Ap[k+1];
            for (; Ai[s] < k && s < s_end; ++s)
                ;
            
            if (s == s_end) // only happens if diagonal is 0, output is not correct in this case.
                A_diag[k] = s-1;
            else
                A_diag[k] = s;
        }
        A_diag_ = A_diag;
        
        return out;
    }
    
    ~CholObj()
    {
        cs_spfree(Lt);
        cs_free(A_diag_);
        cs_free(x_);
        cs_free(c_);
    }
    
    void operator()(const cs<Tv, Ti>* A, cs<Tv, Ti>* out)
    {
        Ti *Ai = A->i, *Ap = A->p; Tv *Ax = A->x;
        Ti nzmax = out->nzmax; Ti n = A->n;
        Ti *Li = out->i, *Lp = out->p;
        Ti *Lti = Lt->i, *Ltp = Lt->p;
        
        const Tv Tv0 = Tv(0); const Tv TvLarge = Tv(1e128);
        Tv *Lx = out->x;
        Tv *x = x_; Ti *c = c_;
        Ti *A_diag = A_diag_;
        
        Ti* Lti_ptr = Lti;
        for (Ti k = 0; k < n; ++k)
        {
            x[k] = Tv0; c[k] = Lp[k];
            
            Ti s_end = A_diag[k];
            for (Ti s = Ap[k]; s < s_end; ++s)
                x[Ai[s]] = Ax[s];
            
            // Solve L_11 l_12 = a_12
            Tv d = Ax[s_end]; Ti i;
            for (; (i = *(Lti_ptr++)) < k;)
            {
                Ti dLi = Lp[i], ci = c[i]++;
                Tv Lki = x[i] / Lx[dLi];
                x[i] = Tv0; // maintain x = 0 for the (k+1) iteration
                
                for (Ti q = dLi + 1; q < ci; ++q)
                    x[Li[q]] -= Lx[q] * Lki;

                d -= Lki * Lki;
                Lx[ci] = Lki;
            }

            // l_22 = sqrt(a22 - <l12,l12>)
            if (d <= Tv0)
                Lx[c[k]++] = TvLarge;
            else
                Lx[c[k]++] = sqrt(d);
        }
    }
};