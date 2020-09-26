#include <algorithm>
#include "CSparse/cs.h"
#include "CSparse/cs_chol.h"
#include "CSparse/cs_schol.h"
#include "CSparse/cs_ereach.h"


template<typename Tv, typename Ti>
struct WeightedMultiplyObj
{
    Tv *tmp;
    
    cs<Tv, Ti>* prepare(const cs<Tv, Ti>* A, const cs<Tv, Ti>* B)
    {
        cs<Tv, Ti>* out = cs_multiply<Tv,Ti>(A, B, false);
        out->x = cs_malloc<Tv>(out->nzmax);
        tmp = cs_malloc<Tv>(out->m);
        
        for (Ti j = 0; j < out->n; j++)
            std::sort(out->i + out->p[j], out->i + out->p[j + 1]);
        
        return out;
    }
    
    ~WeightedMultiplyObj()
    {
        cs_free(tmp);
    }
    
    void operator()(const cs<Tv, Ti>* A, const Tv *w, const cs<Tv, Ti>* B, cs<Tv, Ti>* out)
    {
        Ti m = out->m, n = out->n;
        Ti *Ap = A->p, *Bp = B->p, *Cp = out->p;
        Ti *Ai = A->i, *Bi = B->i, *Ci = out->i;
        Tv *Ax = A->x, *Bx = B->x, *Cx = out->x;
        Tv *x = tmp;
        
        const Tv Tv0 = Tv(0);
        for (Ti j1 = 0; j1 < n; j1++)
        {
            for (Ti p1 = Cp[j1]; p1 < Cp[j1 + 1]; p1++)
                x[Ci[p1]] = Tv0;

            for (Ti p1 = Bp[j1]; p1 < Bp[j1 + 1]; p1++)
            {
                Ti j2 = Bi[p1];
                Tv beta = Bx[p1] * Tv(w[j2]);
                
                for (Ti p2 = Ap[j2]; p2 < Ap[j2 + 1]; p2++)
                {
                    x[Ai[p2]] += beta * Ax[p2];
                }
            }

            for (Ti p1 = Cp[j1]; p1 < Cp[j1 + 1]; p1++)
                Cx[p1] = x[Ci[p1]];
        }
    }
};