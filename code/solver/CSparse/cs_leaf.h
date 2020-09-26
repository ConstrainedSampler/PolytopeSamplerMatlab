// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_LEAF_H_
#define CSPARSE_CS_LEAF_H_
#include "cs.h"

/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) */
namespace CSparse {
	template <typename Ti>
	Ti cs_leaf(Ti i, Ti j, const Ti* first, Ti* maxfirst, Ti* prevleaf,
		Ti* ancestor, Ti* jleaf)
	{
		Ti q, s, sparent, jprev;
		if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1);
		*jleaf = 0;
		if (i <= j || first[j] <= maxfirst[i]) return (-1);  /* j not a leaf */
		maxfirst[i] = first[j];      /* update max first[j] seen so far */
		jprev = prevleaf[i];          /* jprev = previous leaf of ith subtree */
		prevleaf[i] = j;
		*jleaf = (jprev == -1) ? 1 : 2; /* j is first or subsequent leaf */
		if (*jleaf == 1) return (i);   /* if 1st leaf, q = root of ith subtree */
		for (q = jprev; q != ancestor[q]; q = ancestor[q]);
		for (s = jprev; s != q; s = sparent)
		{
			sparent = ancestor[s];    /* path compression */
			ancestor[s] = q;
		}
		return (q);                    /* q = least common ancester (jprev,j) */
	}
}
#endif