#ifndef SIMDVEC_H_
#define SIMDVEC_H_

#include <immintrin.h>

#define REPEAT_1 DECL(0)
#define REPEAT_2 DECL(0) DECL(1)
#define REPEAT_3 DECL(0) DECL(1) DECL(2)
#define REPEAT_4 DECL(0) DECL(1) DECL(2) DECL(3)
#define REPEAT_5 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4)
#define REPEAT_6 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5)
#define REPEAT_7 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6)
#define REPEAT_8 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7)
#define REPEAT_12 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11)
#define REPEAT_16 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15)
#define REPEAT_20 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19)
#define REPEAT_24 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23)
#define REPEAT_28 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23) DECL(24) DECL(25) DECL(26) DECL(27)
#define REPEAT_32 DECL(0) DECL(1) DECL(2) DECL(3) DECL(4) DECL(5) DECL(6) DECL(7) DECL(8) DECL(9) DECL(10) DECL(11) DECL(12) DECL(13) DECL(14) DECL(15) DECL(16) DECL(17) DECL(18) DECL(19) DECL(20) DECL(21) DECL(22) DECL(23) DECL(24) DECL(25) DECL(26) DECL(27) DECL(28) DECL(29) DECL(30) DECL(31)

#define PROGRAM_REPEAT DECL(1) DECL(2) DECL(4) DECL(8) DECL(12) DECL(16) DECL(20) DECL(24) DECL(28) DECL(32)
#define SUPPORTED_MSG "We only support 1, 2, 4, 8, 12, 16, 20, 24, 28 or 32 many vectors."

size_t nearSupportLength(size_t k)
{
    size_t ret = k;
    if (k <= 1)
        return 1;
    
    if (k == 2)
        return 2;
    
    if (k > 32)
        return 32;
    
    return ((k + 3) / 4) * 4;
}

template<typename Tv, size_t k>
struct SIMDVec {};

#ifdef SIMD_DOUBLE
template<typename Tv>
struct ConstVec
{
};

template<>
struct ConstVec<double>
{
    __m256d x;
    __m128d y;
    double z;
    
    ConstVec<double>(double rhs)
    {
        x = _mm256_set1_pd(rhs);
        y = _mm_set1_pd(rhs);
        z = rhs;
    }
};

template<>
struct SIMDVec<double, 1>
{
    double z;
    
    SIMDVec& operator*=(const ConstVec<double> rhs)
    {
        z *= rhs.z;
        return *this;
    }
    
    SIMDVec& sub_assign(const ConstVec<double> a, const SIMDVec<double,1> &b)
    {
        z -= a.z * b.z;
        return *this;
    }
};

template<>
struct SIMDVec<double, 2>
{
    __m128d y;
    
    SIMDVec& operator*=(const ConstVec<double> rhs)
    {
        y = _mm_mul_pd(y, rhs.y);
        
        return *this;
    }
    
    SIMDVec& sub_assign(const ConstVec<double> a, const SIMDVec<double,2> &b)
    {
        y = _mm_fnmadd_pd(a.y, b.y, y);
        
        return *this;
    }
};

#define VEC_K 1
#define VEC_REPEAT REPEAT_1
#include "SIMDVec_double.h"

#define VEC_K 2
#define VEC_REPEAT REPEAT_2
#include "SIMDVec_double.h"

#define VEC_K 3
#define VEC_REPEAT REPEAT_3
#include "SIMDVec_double.h"

#define VEC_K 4
#define VEC_REPEAT REPEAT_4
#include "SIMDVec_double.h"

#define VEC_K 5
#define VEC_REPEAT REPEAT_5
#include "SIMDVec_double.h"

#define VEC_K 6
#define VEC_REPEAT REPEAT_6
#include "SIMDVec_double.h"

#define VEC_K 7
#define VEC_REPEAT REPEAT_7
#include "SIMDVec_double.h"

#define VEC_K 8
#define VEC_REPEAT REPEAT_8
#include "SIMDVec_double.h"

#endif

#ifdef SIMD_TYPE
#define VEC_K 1
#define VEC_REPEAT REPEAT_1
#include "SIMDVec_T.h"

#define VEC_K 2
#define VEC_REPEAT REPEAT_2
#include "SIMDVec_T.h"

#define VEC_K 4
#define VEC_REPEAT REPEAT_4
#include "SIMDVec_T.h"

#define VEC_K 8
#define VEC_REPEAT REPEAT_8
#include "SIMDVec_T.h"

#define VEC_K 12
#define VEC_REPEAT REPEAT_12
#include "SIMDVec_T.h"

#define VEC_K 16
#define VEC_REPEAT REPEAT_16
#include "SIMDVec_T.h"

#define VEC_K 20
#define VEC_REPEAT REPEAT_20
#include "SIMDVec_T.h"

#define VEC_K 24
#define VEC_REPEAT REPEAT_24
#include "SIMDVec_T.h"

#define VEC_K 28
#define VEC_REPEAT REPEAT_28
#include "SIMDVec_T.h"

#define VEC_K 32
#define VEC_REPEAT REPEAT_32
#include "SIMDVec_T.h"

#endif

#endif