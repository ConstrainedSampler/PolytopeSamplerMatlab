#ifndef CSOLVER_H_
#define CSOLVER_H_

#include <random>
#include "CSparse/cs_transpose.h"
#include "CSparse/cs_gaxpy.h"
#include "weighted_multiply.h"
#include "chol.h"
#include "projected_inverse.h"
#include "triple_products.h"
#include "cs_SIMD.h"

namespace env = MexEnvironment;

template <typename Tv, typename Ti>        
struct CSolver
{	// parameters
    cs<Tv,Ti> *A;
    cs<Tv,Ti> *At;
    Tv *w_;
    
    WeightedMultiplyObj<Tv,Ti> weighted_multiply;
    CholObj<Tv,Ti> chol;
    ProjectedInverseObj<Tv,Ti> projected_inverse;
    TripleProductsObj<Tv,Ti> triple_product;
    
    cs<Tv,Ti> *H; // output of weighted_multiply
    cs<Tv,Ti> *L; // output of chol
    cs<Tv,Ti> *S; // output of projected_inverse
    Tv *tau_; // output of triple_product
    Tv* d_; // the JL direction
    Tv *Atd_; // A' * d
    std::mt19937_64 gen;
    
    template <typename Tv_, typename Ti_>
    CSolver(const cs<Tv_, Ti_>* A_, uint64_t seed) : gen(seed)
    {
        A = cs_copy<Tv_, Ti_, Tv, Ti>(A_); // Copy A_ to A
        At = cs_transpose<Tv,Ti>(A, true);
        w_ = cs_malloc<Tv>(A->n);
        H = weighted_multiply.prepare(A, At);
        L = chol.prepare(H);
        S = projected_inverse.prepare(L);
        tau_ = triple_product.prepare(A, S, At);
        d_ = cs_calloc<Tv>(A->m * nearSupportLength(1024));
        Atd_ = cs_calloc<Tv>(A->n * nearSupportLength(1024));
    }
    
    ~CSolver()
    {
        cs_spfree(A);
        cs_spfree(At);
        cs_free(w_);
        cs_spfree(H);
        cs_spfree(L);
        cs_spfree(S);
        cs_free(tau_);
        cs_free(d_);
        cs_free(Atd_);
    }
    
    template <typename Tv_>
    void decompose(const Tv_ *w_in)
    {   // record w
        Ti n = A->n; Tv *w = w_;
        for (Ti j = 0; j < n; j++)
            w[j] = w_in[j];
        
        // compute chol
        weighted_multiply(A, w, At, H);
        chol(H, L);
    }
    
    // nSketch = 0 means we compute the leverage score exactly
    void leverageScoreComplement(size_t nSketch, Tv *out)
    {
        Tv Tv0 = Tv(0.0), Tv1 = Tv(1.0), tau_scale;
        Tv *tau = tau_, *w = w_;
        Ti n = A->n, m = A->m; 
        
        if (nSketch == 0)
        {
            projected_inverse(L, S);

            Tv* Sv = S->x; Ti* Sp = S->p;
            for (Ti k = 0; k < m; ++k)
                Sv[Sp[k]] = Sv[Sp[k]] / Tv(2.0);   // divide the diagonal by 2

            triple_product(A, S, At, tau_);
            tau_scale = Tv(2.0);
        }
        else
        {
            size_t nDim = 0;
            int sketchLeft = (int)nSketch;
            while (sketchLeft > 0)
            {
                size_t k = nearSupportLength(sketchLeft);
                
                #define DECL(k) case k: leverageScoreJL<k>(tau); break;
                switch(k)
                {
                PROGRAM_REPEAT;
                default:
                    env::error(SUPPORTED_MSG);
                }
                #undef DECL
                
                sketchLeft -= int(k);
                nDim += k;
            }
            tau_scale = Tv(1.0/nDim);
        }
        
        for (Ti j = 0; j < n; j++)
        {
            out[j] = Tv1 - tau_scale * tau[j] * w[j];
            tau[j] = Tv0;
        }
    }
    
    template<size_t k>
    void leverageScoreJL(Tv *out_)
    {
        Ti m = A->m, n = A->n;
        Tv *d = randSign(m * k), *Atd = Atd_, *out = out_;
        Lt_solve<Tv, Ti, k>(L, d, d);
        cs_gaxpy<Tv, Ti, k>(At, d, Atd);
        
        Tv Tv0 = Tv(0.0);
        for (Ti i = 0; i < n; ++i)
        {
            Tv ret_i = Tv0;
            for (Ti j = 0; j < k; ++j)
            {
                ret_i += Atd[i*k+j] * Atd[i*k+j];
                Atd[i*k+j] = Tv0;
            }
            
            out[i] += ret_i;
        }
    }
    
    Tv logdet()
    {
        Ti *Li = L->i, *Lp = L->p; Tv *Lx = L->x;
        Ti m = A->m;
        
        Tv ret = Tv(0);
        for (Ti j = 0; j < m; j++)
            ret += log(Lx[Lp[j]]);
        
        return ret * Tv(2.0);
    }

    void diagL(Tv* out)
    {
        Ti* Li = L->i, * Lp = L->p; Tv* Lx = L->x;
        Ti m = A->m;

        for (Ti j = 0; j < m; j++)
            out[j] = Lx[Lp[j]];
    }
    
    template <size_t k>
    void solve(const Tv *b, Tv *out)
    {
        L_solve<Tv, Ti, k>(L, b, out);
        Lt_solve<Tv, Ti, k>(L, out, out);
    };
    
    // This is for internal use
    Tv *randSign(Ti n)
    {
        Tv *d = d_;
        
        for (Ti i = 0; i < n; ++i)
            d[i] = Tv(2.0 * (gen() & 1) - 1.0);
        
        return d;
    }
};



int main()
{
    const char* cmd = env::inputString();
    uint64_t uid = env::inputScalar<uint64_t>();
    if (!strcmp(cmd, "init"))
    {
        auto A = env::inputSparseArray<double>();
        auto solver = new CSolver<Scalar, Index>(&A, uid);
        env::outputScalar<uint64_t>((uint64_t)solver);
    }
    else
    {
        CSolver<Scalar, env::mexIdx>* solver = (CSolver<Scalar, env::mexIdx> *)uid;
        size_t n = solver->A->n; size_t m = solver->A->m;
        if (!strcmp(cmd, "solve"))
        {
            size_t k = -1;
            Scalar *b; Scalar *out;
            double *b_double = (double*)env::inputArray<double>(k, m);
            if (std::is_same<Scalar, double>::value)
            {
                b = (Scalar *)b_double;
                out = (Scalar *)env::outputArray<double>(k, m);
            }
            else
            {
                b = new Scalar[k * m];
                for (size_t i = 0; i < k * m; ++i)
                    b[i] = convert<double, Scalar>(b_double[i]);
                
                out = new Scalar[k * m];
            }
            
            #define DECL(k) case k: solver->solve<k>(b, out); break;
            switch(k)
            {
                PROGRAM_REPEAT;
                default:
                    env::error(SUPPORTED_MSG);
            }
            #undef DECL
            
            if (!std::is_same<Scalar, double>::value)
            {
                env::outputDoubleArray<Scalar>(out, k, m);
                delete[] b;
                delete[] out;
            }
        }
        else if (!strcmp(cmd, "decompose"))
        {
            auto w = env::inputArray<double>(n);
            solver->decompose(w);
        }
        else if (!strcmp(cmd, "leverageScoreComplement"))
        {
            auto k = (size_t)env::inputScalar<double>();
            
            if (std::is_same<Scalar, double>::value)
            {
                auto out = env::outputArray<double>(n);
                solver->leverageScoreComplement(k, (Scalar*)out);
            }
            else
            {
                Scalar *out = new Scalar[n];
                solver->leverageScoreComplement(k, out);
                env::outputDoubleArray<Scalar>(out, n);
                delete[] out;
            }
        }
        else if (!strcmp(cmd, "logdet"))
        {
            env::outputScalar<double>(convert<Scalar, double>(solver->logdet()));
        }
        else if (!strcmp(cmd, "diagL"))
        {
            if (std::is_same<Scalar, double>::value)
            {
                auto out = env::outputArray<double>(m);
                solver->diagL((Scalar*)out);
            }
            else
            {
                Scalar* out = new Scalar[m];
                solver->diagL(out);
                env::outputDoubleArray<Scalar>(out, m);
                delete[] out;
            }
        }
        else if (!strcmp(cmd, "delete"))
        {
            delete solver;
        }
        else
            env::error("Invalid operation.");
    }
}
#endif