#include "qd/dd_real.h"
#include "mex_utils.h"

#define SIMD_TYPE
typedef dd_real Scalar;
typedef MexEnvironment::mexIdx Index;

template<>
double convert(dd_real x)
{
    return to_double(x);
}

#include "CSolver.h"