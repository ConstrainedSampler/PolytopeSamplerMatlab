#include <immintrin.h>
#include "SIMDVec.h"

// x and out can be the same pointer
template<typename Tv, typename Ti, size_t k>
void L_solve(const cs<Tv, Ti>* L, const Tv* x, Tv* out_)
{
    Tv* Lx = L->x;
    Ti m = L->m, n = L->n, *Lp = L->p, *Li = L->i;
    if (x != out_)
        std::copy(x, x + n * k, out_);

    SIMDVec<Tv, k> *out = (SIMDVec<Tv, k> *)out_;

    for (Ti j = 0; j < n; j++)
    {
        //out[j] /= Lx[Lp[j]];
        SIMDVec<Tv, k> out_j = out[j];
        out_j *= 1/Lx[Lp[j]];
        out[j] = out_j;
        
        for (Ti p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {   //out[Li[p]] -= Lx[p] * out[j];
            out[Li[p]].sub_assign(Lx[p], out_j);
        }
    }
}

// x and out can be the same pointer
template<typename Tv, typename Ti, size_t k>
void Lt_solve(const cs<Tv, Ti>* L, const Tv* x, Tv* out_)
{
    Tv* Lx = L->x;
    Ti m = L->m, n = L->n, *Lp = L->p, *Li = L->i;
    if (x != out_)
        std::copy(x, x + m * k, out_);
    
    SIMDVec<Tv, k> *out = (SIMDVec<Tv, k> *)out_;
    for (Ti j = n - 1; j >= 0; j--)
    {
        SIMDVec<Tv, k> out_j = out[j];

        for (Ti p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {   //out[j] -= Lx[p] * out[Li[p]];
            out_j.sub_assign(Lx[p], out[Li[p]]);
        }

        //out_[j] = oj / Lx[Lp[j]];
        out_j *= 1/Lx[Lp[j]];
        out[j] = out_j;
    }
}

template<typename Tv, typename Ti, size_t k>
void cs_gaxpy(const cs<Tv, Ti>* A, const Tv* x_, Tv* y_)
{
    Tv* Ax = A->x;
    Ti m = A->m, n = A->n, *Ap = A->p, *Ai = A->i;
    SIMDVec<Tv, k> *y = (SIMDVec<Tv, k> *)y_;
    SIMDVec<Tv, k> *x = (SIMDVec<Tv, k> *)x_;
    
    for (Ti j = 0; j < n; j++)
    {
        SIMDVec<Tv, k> xj = x[j];
        for (Ti p = Ap[j]; p < Ap[j + 1]; p++)
        {   //y[Ai[p]] += Ax[p] * x[j];
            y[Ai[p]].sub_assign(Ax[p], xj);
        }
    }
}