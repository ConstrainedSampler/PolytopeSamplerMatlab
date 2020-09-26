#include "CSparse/cs_chol.h"
#include "CSparse/cs_schol.h"
#include "CSparse/cs_ereach.h"

// Compute inv(L L') restricted on L
template<typename Tv, typename Ti>
struct ProjectedInverseObj
{
    cs<bool, Ti> *Lt = nullptr;
    
    // workspaces
    Tv *x_ = nullptr; 
    Ti *c_ = nullptr; // c[i] = index the last nonzero on column i in the current L
    
    cs<Tv, Ti>* prepare(const cs<Tv, Ti>* L)
    {
        Ti n = L->n;
        c_ = cs_malloc<Ti>(n);
        x_ = cs_malloc<Tv>(n);
        
        cs<Tv, Ti>* out = cs_copy<Tv, Ti, Tv, Ti>(L);
        Lt = cs_transpose<bool,Ti>((cs<bool, Ti>*)L, false);
        
        return out;
    }
    
    ~ProjectedInverseObj()
    {
        cs_spfree(Lt);
        cs_free(x_);
        cs_free(c_);
    }
    
    void operator()(const cs<Tv, Ti>* L, cs<Tv, Ti>* out)
    {
        Tv* Sv = out->x; Ti n = out->n;
        Ti* Li = L->i, *Lp = L->p; Tv* Lv = L->x;
        Ti* Lti = Lt->i, *Ltp = Lt->p;
        Tv* x = x_;
        Ti* c = c_;
        Tv Tv0 = Tv(0);
        
        for (Ti k = 0; k < n; ++k)
            c[k] = Lp[k + 1] - 1;

        for (Ti k = n - 1; k != -1; --k)
        {
            for (Ti p = Lp[k] + 1; p < Lp[k+1]; ++p)
                x[Li[p]] = Sv[p];

            Tv sum = 1 / Lv[Lp[k]];
            for (Ti p = Ltp[k + 1] - 1; p != Ltp[k] - 1; --p)
            {
                Ti i = Lti[p], Lpi = Lp[i];

                for (Ti q = Lp[i + 1] - 1; q > Lpi; q--)
                    sum -= Lv[q] * x[Li[q]];
                
                sum = sum / Lv[Lpi];
                x[i] = sum;
                Sv[c[i]] = sum;
                c[i]--;
                sum = Tv0;
            }
        }
    }
};

