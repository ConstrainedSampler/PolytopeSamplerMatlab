template<>
struct SIMDVec<double, VEC_K * 4>
{
    //__m256d x[VEC_K];
#define DECL(k) __m256d x ## k;
        VEC_REPEAT;
#undef DECL
    
    SIMDVec& operator*=(const ConstVec<double> rhs)
    {
        //for (int i = 0; i < VEC_K; ++i)
        //    x[i] = _mm256_mul_pd(x[i], rhs.x);
        
#define DECL(k) x ## k = _mm256_mul_pd(x ## k, rhs.x);
        VEC_REPEAT;
#undef DECL
        return *this;
    }
    
    SIMDVec& sub_assign(const ConstVec<double> a, const SIMDVec<double,VEC_K * 4> &b)
    {
        //for (int i = 0; i < VEC_K; ++i)
        //    x[i] = _mm256_fnmadd_pd(a.x, b.x[i], x[i]);
        
#define DECL(k) x ## k = _mm256_fnmadd_pd(a.x, b.x ## k, x ## k);
        VEC_REPEAT;
#undef DECL
        return *this;
    }
};

#undef VEC_K
#undef VEC_REPEAT
