// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_DROPTOL_H_
#define CSPARSE_CS_DROPTOL_H_
#include "cs.h"

namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_tol(Ti i, Ti j, Tv aij, void* tol)
	{
		return (fabs(aij) > * ((Tv*)tol));
	}

	template <typename Tv, typename Ti>
	Ti cs_droptol(cs<Tv,Ti>* A, Tv tol)
	{
		return (cs_fkeep(A, &cs_tol, &tol));    /* keep all large entries */
	}
}

#endif