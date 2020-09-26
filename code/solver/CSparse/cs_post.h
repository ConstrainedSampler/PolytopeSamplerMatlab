// Modified from CSparse 3.2. License: LGPL-3.0
#ifndef CSPARSE_CS_POST_H_
#define CSPARSE_CS_POST_H_
#include "cs.h"
#include "cs_tdfs.h"

/* post order a forest */
namespace CSparse {
	template <typename Ti>
	Ti* cs_post(const Ti* parent, Ti n)
	{
		Ti j, k = 0, * post, * w, * head, * next, * stack;
		if (!parent) return (NULL);                        /* check inputs */
		post = cs_malloc<Ti>(n);                /* allocate result */
		w = cs_malloc<Ti>(3 * n);                 /* get workspace */
		if (!w || !post) return (cs_idone<Ti,Ti>(post, NULL, w, 0));
		head = w; next = w + n; stack = w + 2 * n;
		for (j = 0; j < n; j++) head[j] = -1;           /* empty linked lists */
		for (j = n - 1; j >= 0; j--)            /* traverse nodes in reverse order*/
		{
			if (parent[j] == -1) continue;    /* j is a root */
			next[j] = head[parent[j]];      /* add j to list of its parent */
			head[parent[j]] = j;
		}
		for (j = 0; j < n; j++)
		{
			if (parent[j] != -1) continue;    /* skip j if it is not a root */
			k = cs_tdfs(j, k, head, next, post, stack);
		}
		return (cs_idone<Ti, Ti>(post, NULL, w, 1));  /* success; free w, return post */
	}
}
#endif
