#include <unordered_map>
#include "mex.h"
#include "cpu_x86.h"
#include "cpu_x86.cpp"

using namespace std;
using namespace FeatureDetector;

template <typename T>
mxArray* mxCreateScalar(T x);

template <>
mxArray* mxCreateScalar(bool x)
{
	return mxCreateLogicalScalar(x);
}

template <typename T>
mxArray* createStructure(const unordered_map<const char*, T> &map)
{
	const char** fieldnames = new const char* [map.size()];
	mxArray** values = new mxArray *[map.size()];

	size_t k = 0;

	for (const auto& nv : map)
	{
		fieldnames[k] = (char*)mxMalloc(64);
		values[k] = mxCreateScalar<T>(nv.second);
		memcpy((void*)fieldnames[k], nv.first, strlen(nv.first) + 1);
		++k;
	}

	auto pt = mxCreateStructMatrix(1, 1, (int)map.size(), fieldnames);
	for (size_t i = 0; i < map.size(); ++i)
	{
		mxFree((void*)fieldnames[i]);
		mxSetFieldByNumber(pt, 0, i, values[i]);
	}

	delete[] fieldnames;
	delete[] values;

	return pt;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nlhs == 1)
	{
		cpu_x86 features;
		features.detect_host();

		unordered_map<const char*, bool> info;
		info["SSE"] = features.HW_SSE;
		info["SSE2"] = features.HW_SSE2;
		info["SSE3"] = features.HW_SSE3;
		info["SSE41"] = features.HW_SSE41;
		info["SSE42"] = features.HW_SSE42;
		info["AVX"] = features.HW_AVX;
		info["AVX2"] = features.HW_AVX2;
		info["FMA"] = features.HW_FMA3;
		info["AVX512F"] = features.HW_AVX512_F;
		info["OS_AVX"] = features.OS_AVX;
		info["OS_AVX512"] = features.OS_AVX512;
      
		plhs[0] = createStructure(info);
	}
}
