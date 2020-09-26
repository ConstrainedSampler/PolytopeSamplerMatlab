// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_PINV_H_
#define CSPARSE_CS_PINV_H_
#include "cs.h"

/* pinv = p', or p = pinv' */
namespace CSparse {
	template <typename Ti>
	Ti* cs_pinv(Ti const* p, Ti n)
	{
		Ti k, * pinv;
		if (!p) return (NULL);                     /* p = NULL denotes identity */
		pinv = cs_malloc<Ti>(n);        /* allocate result */
		if (!pinv) return (NULL);                  /* out of memory */
		for (k = 0; k < n; k++) pinv[p[k]] = k;/* invert the permutation */
		return (pinv);                             /* return result */
	}
}

#endif