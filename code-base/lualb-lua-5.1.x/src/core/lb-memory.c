#include <stdlib.h>
#include <string.h>
#include "lb-memory.h"
#include "lb-messages.h"

void*
lb_malloc(size_t size)
{
	void* ptr = malloc(size);

	if (ptr == NULL)
		lb_fatal("failed to allocate %d bytes of memory\n", size);

	return ptr;
}

void
lb_free(void* ptr)
{
	if (ptr)
		free(ptr);
}

void*
lb_realloc(void* ptr, size_t size)
{
	return realloc(ptr, size);
}

void
lb_scpy(double*       LB_RESTRICT dst,
	  const double* LB_RESTRICT src,
	  int rank, int stride, int count)
{
	lb_assert(dst && src);
	lb_assert(rank > 0 && stride != 0);

	rank *= sizeof(double);
	for (; count > 0; --count) {
		memcpy(dst, src, rank);
		dst += stride;
		src += stride;
	}
}
