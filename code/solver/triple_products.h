// Compute diag(A S B)
template<typename Tv, typename Ti>
struct TripleProductsObj
{
    Tv *out;
    Tv *s_col_;
    Ti *s_mark_;
    
    Tv *prepare(const cs<Tv, Ti>* A, const cs<Tv, Ti>* S, const cs<Tv, Ti>* B)
    {
        s_col_ = cs_malloc<Tv>(A->m);
        s_mark_ = cs_malloc<Ti>(A->m);
        return cs_calloc<Tv>(A->n);
    }
    
    ~TripleProductsObj()
    {
        cs_free(s_col_);
        cs_free(s_mark_);
    }
    
    void operator()(const cs<Tv, Ti>* A, const cs<Tv, Ti>* S, const cs<Tv, Ti>* B, Tv *out_)
    {
        Ti m = A->m, n = B->n, k = A->n;
        Ti *Ai = A->i, *Ap = A->p, *Bi = B->i, *Bp = B->p, *Si = S->i, *Sp = S->p;
        Tv *Ax = A->x, *Bx = B->x, *Sx = S->x;
        Tv *s_col = s_col_;
        Ti *s_mark = s_mark_;
        Tv *out = out_;
        
        std::fill(s_mark, s_mark + m, Ti(-1));

        for (Ti j = 0; j < n; j++)
        {
            for (Ti p = Sp[j]; p < Sp[j + 1]; p++)
            {
                s_col[Si[p]] = Sx[p];
                s_mark[Si[p]] = j;
            }

            for (Ti p = Bp[j]; p < Bp[j + 1]; p++)
            {
                Ti i = Bi[p]; Tv b = Bx[p];
                for (Ti q = Ap[i]; q < Ap[i + 1]; q++)
                {
                    Tv a = Ax[q]; Ti a_i = Ai[q];
                    if (s_mark[a_i] == j)
                        out[i] += s_col[a_i] * a * b;
                }
            }
        }
    }

};
