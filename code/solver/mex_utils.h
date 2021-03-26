#pragma once
#include "mex.h"
#include <cstdint>
#include <typeinfo>
#include <string>
#include "PackedCSparse/SparseMatrix.h"
using namespace PackedCSparse;

namespace MexEnvironment
{
    typedef std::make_signed<mwIndex>::type mexIdx;
    
    template <typename T> mxClassID mexType()
    {
        return mxUNKNOWN_CLASS;
    }
    
    #define MEX_TYPE_DEFINE(CType, MexType) \
        template <> mxClassID mexType<CType>() { return MexType; }
    
    MEX_TYPE_DEFINE(bool, mxLOGICAL_CLASS);
    MEX_TYPE_DEFINE(char, mxCHAR_CLASS);
    MEX_TYPE_DEFINE(int8_t, mxINT8_CLASS);
    MEX_TYPE_DEFINE(uint8_t, mxUINT8_CLASS);
    MEX_TYPE_DEFINE(int16_t, mxINT16_CLASS);
    MEX_TYPE_DEFINE(uint16_t, mxUINT16_CLASS);
    MEX_TYPE_DEFINE(int32_t, mxINT32_CLASS);
    MEX_TYPE_DEFINE(uint32_t, mxUINT32_CLASS);
    MEX_TYPE_DEFINE(int64_t, mxINT64_CLASS);
    MEX_TYPE_DEFINE(uint64_t, mxUINT64_CLASS);
    MEX_TYPE_DEFINE(float, mxSINGLE_CLASS);
    MEX_TYPE_DEFINE(double, mxDOUBLE_CLASS);
	
	#undef MEX_TYPE_DEFINE
    
    /* ====== Parameters Info ====== */
    // Input
    const mxArray **prhs;
    size_t nrhs;
    size_t rhs_id; // pointer to next input
    
    // Output
    mxArray **plhs;
    size_t nlhs;
    size_t lhs_id; // pointer to next output
    
    const mxArray *input()
    {
        if (rhs_id >= nrhs)
            throw "At least " + std::to_string(rhs_id + 1) + " input parameters are required.";
        
        return prhs[rhs_id++];
    }
    
    void output(mxArray *out)
    {
        if (lhs_id >= nlhs)
            throw "At least " + std::to_string(lhs_id + 1) + " output parameters are required.";
        
        plhs[lhs_id++] = out;
    }
    
    void errorSize(size_t mRequired, size_t nRequired, size_t m, size_t n)
    {
        if (mRequired == -1 || nRequired == -1)
            throw "Incorrect dimension for " + std::to_string(rhs_id) + "-th parameter. "
                    "It should be (" + std::to_string(mRequired) + "," + std::to_string(nRequired) + ")"
                    " where -1 indicates any non-negative numbers.";
        else
            throw "Incorrect dimension for " + std::to_string(rhs_id) + "-th parameter. "
                    "It should be (" + std::to_string(mRequired) + "," + std::to_string(nRequired) + ")"
                    " instead of (" + std::to_string(m) + "," + std::to_string(n) + ")";
    }
    
    /* ====== Various Input Functions ====== */
    template<typename T>
    const T *inputArray(size_t &mRequired, size_t &nRequired, const mxArray *pt = nullptr)
    {
        if (pt == nullptr)
            pt = input();
        
        auto nDim = mxGetNumberOfDimensions(pt);
        
        if (mxIsComplex(pt) || mxIsSparse(pt))
            throw "The " + std::to_string(rhs_id) + "-th parameter should be a real full array.";
        
        auto dims = mxGetDimensions(pt);
        if (mxGetClassID(pt) != mexType<T>())
            throw "The " + std::to_string(rhs_id) + "-th parameter should be " + typeid(T).name();
        
        size_t m = dims[0], n = dims[1];
        
        if ((mRequired != -1 && mRequired != m) || (nRequired != -1 && nRequired != n))
            errorSize(mRequired, nRequired, m, n);
        mRequired = m; nRequired = n;
        
        return (T*)mxGetData(pt);
    }
    
    template<typename T>
    const T *inputArray(size_t &m)
    {
        size_t n = 1;
        return inputArray<T>(m, n);
    }
    
    template<typename T>
    T inputScalar()
    {
        size_t m = 1, n = 1;
        return *inputArray<T>(m, n);
    }

    template<typename T>
    T inputScalar(T default_value)
    {
        T val = default_value;
        if (rhs_id < nrhs)
            val = inputScalar<T>();
        return val;
    }
    
    const char *inputString()
    {
        const mxArray *pt = input();
        if (mxIsChar(pt) != 1 || mxGetM(pt)!=1)
            throw "The " + std::to_string(rhs_id) + "-th parameter should be a string.";
        
        return mxArrayToString(pt);
    }

    //A better thing to do is to have a matrix "view" and output the view
    template<typename Tv = double, typename Ti = mexIdx>
    SparseMatrix<Tv, Ti> inputSparseArray(size_t mRequired = -1, size_t nRequired = -1)
    {
        const mxArray* pt = input();
        size_t m, n, nzmax;
        Tv* x; mexIdx* ir, * jc;
        if (mxIsComplex(pt) || mxGetNumberOfDimensions(pt) != 2 || !mxIsSparse(pt))
            throw "The " + std::to_string(rhs_id) + "-th parameter should be a real sparse array.";
        
        if (mxGetClassID(pt) != mexType<Tv>())
            throw "The " + std::to_string(rhs_id) + "-th parameter should be " + typeid(Tv).name();
        
        m = mxGetM(pt); n = mxGetN(pt); nzmax = mxGetNzmax(pt);
        x = (Tv*)mxGetData(pt);
        ir = (mexIdx*)mxGetIr(pt);
        jc = (mexIdx*)mxGetJc(pt);
            
        if ((mRequired != -1 && mRequired != m) || (nRequired != -1 && nRequired != n))
            errorSize(mRequired, nRequired, m, n);
        
        if (!(std::is_same<Tv, double>::value || std::is_same<Tv, bool>::value))
            throw "The " + std::to_string(rhs_id) + "-th parameter should be " + typeid(Tv).name();
        
        
        SparseMatrix<Tv, Ti> A(m, n, nzmax);
        std::copy(ir, ir + nzmax, A.i.get());
        std::copy(jc, jc + n + 1, A.p.get());
        std::copy(x, x + nzmax, A.x.get());
        return std::move(A);
    }
    
    /* ====== Various Output Functions ====== */
    template<typename Tv = double, typename Ti = mexIdx>
    void outputSparseArray(const SparseMatrix<Tv, Ti>& A, bool patternOnly = false)
    {
        mxArray* pt;
        Ti nz = A.nnz();
		if (patternOnly)
		{
			pt = mxCreateSparseLogicalMatrix(A.m, A.n, nz);
			bool* Ax = (bool*)mxGetData(pt);
			for (Ti s = 0; s < nz; ++s)
				Ax[s] = bool(A.x[s]);
		}
		else
		{
			pt = mxCreateSparse(A.m, A.n, nz, mxREAL);
			double* Ax = (double*)mxGetData(pt);
			for (Ti s = 0; s < nz; ++s)
				Ax[s] = double(A.x[s]);
		}
        output(pt);

        mexIdx* Ai = (mexIdx*)mxGetIr(pt);
        for (Ti s = 0; s < nz; ++s)
            Ai[s] = mexIdx(A.i[s]);

        mexIdx* Ap = (mexIdx*)mxGetJc(pt);
        for (Ti s = 0; s <= A.n; ++s)
            Ap[s] = mexIdx(A.p[s]);
    }

    template<typename Tv>
    Tv* outputArray(size_t m, size_t n = 1)
    {
		mxArray* pt; size_t l = 1;
		Tv* d = (Tv*)mxMalloc((m * n * l) * sizeof(Tv));
		pt = mxCreateNumericMatrix(0, 0, mexType<Tv>(), mxREAL);

		mxSetData(pt, d);
        output(pt);
        mwSize dims[3] = { (mwSize)m, (mwSize)n, (mwSize)l };
        mxSetDimensions(pt, dims, 3);
        return d;
    }

	template<typename Tv>
	void outputDoubleArray(Tv *x, size_t m, size_t n = 1)
	{
		double *out = outputArray<double>(m, n);
		for (size_t s = 0; s < m * n; ++s)
			out[s] = double(x[s]);
	}
    
    template<typename T>
    void outputScalar(T val)
    {
        *outputArray<T>(1, 1) = val;
    }
};

int main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try
    {
        MexEnvironment::nlhs = nlhs;
        MexEnvironment::nrhs = nrhs;
        MexEnvironment::plhs = plhs;
        MexEnvironment::prhs = prhs;
        MexEnvironment::lhs_id = 0;
        MexEnvironment::rhs_id = 0;
        main();
    }
    catch (const char *str)
    {
        mexErrMsgTxt(str);
    }
    catch (std::string str)
    {
        mexErrMsgTxt(str.c_str());
    }
    catch (const std::exception& e)
    {
        mexErrMsgTxt(e.what());
    }
    return;
}