// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_PVEC_H_
#define CSPARSE_CS_PVEC_H_
#include "cs.h"

/* x = b(p), for dense vectors x and b; p=NULL denotes identity */
namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_pvec(const Ti* p, const Tv* b, Tv* x, Ti n)
	{
		Ti k;
		if (!x || !b) return (0);                              /* check inputs */
		for (k = 0; k < n; k++) x[k] = b[p ? p[k] : k];
		return (1);
	}
}
#endif