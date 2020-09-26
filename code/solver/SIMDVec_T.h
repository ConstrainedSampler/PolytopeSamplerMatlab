template<typename T>
struct SIMDVec<T, VEC_K>
{
    //T x[VEC_K];
#define DECL(k) T x ## k;
        VEC_REPEAT;
#undef DECL
    
    SIMDVec& operator*=(const T rhs)
    {
        //for (int i = 0; i < VEC_K; ++i)
        //    x[i] *= rhs;
        
#define DECL(k) x ## k *= rhs;
        VEC_REPEAT;
#undef DECL
        return *this;
    }
    
    SIMDVec& sub_assign(const T a, const SIMDVec<T, VEC_K> &b)
    {
        //for (int i = 0; i < VEC_K; ++i)
        //    x[i] -= a * b.x[i];
        
#define DECL(k) x ## k -= a * b.x ## k;
        VEC_REPEAT;
#undef DECL
        return *this;
    }
};

#undef VEC_K
#undef VEC_REPEAT
