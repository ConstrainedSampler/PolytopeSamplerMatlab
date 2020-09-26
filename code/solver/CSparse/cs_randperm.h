// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_RANDPERM_H_
#define CSPARSE_CS_RANDPERM_H_
#include "cs.h"

/* return a random permutation vector, the identity perm, or p = n-1:-1:0.
 * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
 * p = random permutation.  */
namespace CSparse {
	template <typename Ti>
	Ti* cs_randperm(Ti n, Ti seed)
	{
		Ti* p, k, j, t;
		if (seed == 0) return (NULL);      /* return p = NULL (identity) */
		p = cs_malloc<Ti>(n);   /* allocate result */
		if (!p) return (NULL);             /* out of memory */
		for (k = 0; k < n; k++) p[k] = n - k - 1;
		if (seed == -1) return (p);        /* return reverse permutation */
		srand((int)seed);                      /* get new random number seed */
		for (k = 0; k < n; k++)
		{
			j = k + (rand() % (n - k));    /* j = rand integer in range k to n-1 */
			t = p[j];                     /* swap p[k] and p[j] */
			p[j] = p[k];
			p[k] = t;
		}
		return (p);
	}
}
#endif