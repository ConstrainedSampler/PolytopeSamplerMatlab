// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_HOUSE_H_
#define CSPARSE_CS_HOUSE_H_
#include "cs.h"

/* create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*e1.  See Algo 5.1.1, Golub & Van Loan, 3rd ed. */
namespace CSparse {
	template <typename Tv, typename Ti>
	Tv cs_house(Tv* x, Tv* beta, Ti n)
	{
		Tv s, sigma = 0;
		Ti i;
		if (!x || !beta) return (-1);          /* check inputs */
		for (i = 1; i < n; i++) sigma += x[i] * x[i];
		if (sigma == 0)
		{
			s = fabs(x[0]);                  /* s = |x(0)| */
			(*beta) = (x[0] <= 0) ? 2 : 0;
			x[0] = 1;
		}
		else
		{
			s = sqrt(x[0] * x[0] + sigma);  /* s = norm (x) */
			x[0] = (x[0] <= 0) ? (x[0] - s) : (-sigma / (x[0] + s));
			(*beta) = -1. / (s * x[0]);
		}
		return (s);
	}
}
#endif
