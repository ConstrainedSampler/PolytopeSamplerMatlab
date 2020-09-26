// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_FKEEP_H_
#define CSPARSE_CS_FKEEP_H_
#include "cs.h"

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
namespace CSparse {
	template <typename Tv, typename Ti>
	Ti cs_fkeep(cs<Tv,Ti>* A, bool(*fkeep) (Ti, Ti, Tv, void*), void* other, bool reallocate = true)
	{
		Ti j, p, nz = 0, n, * Ap, * Ai;
		Tv* Ax;
		if (!CS_CSC(A) || !fkeep) return (-1);    /* check inputs */
		n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
		for (j = 0; j < n; j++)
		{
			p = Ap[j];                        /* get current location of col j */
			Ap[j] = nz;                       /* record new location of col j */
			for (; p < Ap[j + 1]; p++)
			{
				if (fkeep(Ai[p], j, Ax ? Ax[p] : 1, other))
				{
					if (Ax) Ax[nz] = Ax[p];  /* keep A(i,j) */
					Ai[nz++] = Ai[p];
				}
			}
		}
		Ap[n] = nz;                           /* finalize A */
        if (reallocate)
    		cs_sprealloc<Tv,Ti>(A, 0);                   /* remove extra space from A */
		return (nz);
	}
}

#endif