// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_DROPZEROS_H_
#define CSPARSE_CS_DROPZEROS_H_
#include "cs.h"
#include "cs_fkeep.h"

namespace CSparse {
	template <typename Tv, typename Ti>
	bool cs_nonzero(Ti i, Ti j, Tv aij, void* other)
	{
		return (aij != 0);
	}

	template <typename Tv, typename Ti>
	Ti cs_dropzeros(cs<Tv,Ti>* A, bool reallocate = true)
	{
		return (cs_fkeep(A, &cs_nonzero<Tv,Ti>, NULL, reallocate));  /* keep all nonzero entries */
	}
}

#endif