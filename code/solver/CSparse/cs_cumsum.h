// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_CUMSUM_H_
#define CSPARSE_CS_CUMSUM_H_
#include "cs.h"
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */

namespace CSparse {
	template <typename Ti>
	int64_t cs_cumsum(Ti*p, Ti*c, Ti n)
	{
		Ti i, nz = 0 ;
		int64_t nz2 = 0 ;
		if (!p || !c) return (-1) ;     /* check inputs */
		for (i = 0 ; i < n ; i++)
		{
			p [i] = nz ;
			nz += c [i] ;
			nz2 += (int64_t)c [i] ;              /* also in double to avoid csi overflow */
			c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
		}
		p [n] = nz ;
		return (nz2) ;                  /* return sum (c [0..n-1]) */
	}
}

#endif