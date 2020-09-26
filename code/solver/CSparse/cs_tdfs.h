// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_TDFS_H_
#define CSPARSE_CS_TDFS_H_
#include "cs.h"

/* depth-first search and postorder of a tree rooted at node j */
namespace CSparse {
	template <typename Ti>
	Ti cs_tdfs(Ti j, Ti k, Ti* head, const Ti* next, Ti* post, Ti* stack)
	{
		Ti i, p, top = 0;
		if (!head || !next || !post || !stack) return (-1);    /* check inputs */
		stack[0] = j;                 /* place j on the stack */
		while (top >= 0)                /* while (stack is not empty) */
		{
			p = stack[top];           /* p = top of stack */
			i = head[p];              /* i = youngest child of p */
			if (i == -1)
			{
				top--;                 /* p has no unordered children left */
				post[k++] = p;        /* node p is the kth postordered node */
			}
			else
			{
				head[p] = next[i];   /* remove i from children of p */
				stack[++top] = i;     /* start dfs on child node i */
			}
		}
		return (k);
	}
}
#endif
