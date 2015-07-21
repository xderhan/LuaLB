#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "lb-memory.h"
#include "lb-messages.h"
#include "lb-cart-partition.h"

static int
brick_face_area(int n_dims, const int* size, int face)
{
	int d, area;

	lb_assert(n_dims > 0);

	lb_assert(size != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(size[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(face < n_dims);

	for (d = 0, area = 1; d < n_dims; ++d)
		if (d != face)
			area *= size[d];

	return area;
}

/*****************************************************************************/

static int
brick_communication_cost(int        n_dims,
			 const int* dims,
			 const int* gsize,
			 const int* periods)
{
	int  d, area, done;
	int* coords;
	int* brick_size;

	lb_assert(n_dims > 0);

	lb_assert(dims != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(dims[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(gsize != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(gsize[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(periods != NULL);

	coords = lb_newa(int, n_dims);
	brick_size = lb_newa(int, n_dims);

	for (d = 0; d < n_dims; ++d)
		coords[d] = 0;

	for (area = 0, done = 0; done == 0;) {
		for (d = 0; d < n_dims; ++d)
			if (coords[d] != dims[d] - 1) {
				++coords[d];
				break;
			} else {
				coords[d] = 0;
				if (d == n_dims - 1)
					done = 1;
			}

		lb_cart_partition_size(n_dims, dims, gsize,
					 coords, brick_size);

		for (d = 0; d < n_dims; ++d) {
			int face_area = brick_face_area(n_dims, brick_size, d);

			if (coords[d] != 0 || periods[d])
				area += face_area;

			if (coords[d] != dims[d] - 1 || periods[d])
				area += face_area;
		}
	}

	return area;
}

/*****************************************************************************/

void
lb_cart_partition_size(int        n_dims,
			 const int* dims,
			 const int* gsize,
			 const int* coords,
			 int*       lsize)
{
	int d;

	lb_assert(n_dims > 0);

	lb_assert(dims != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(dims[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(gsize != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(gsize[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(coords != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(coords[d] >= 0 && coords[d] < dims[d]);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(lsize != NULL);

	for (d = 0; d < n_dims; ++d) {
		lsize[d] = gsize[d]/dims[d];

		if (coords[d] == dims[d] - 1)
			lsize[d] += gsize[d] % dims[d];
	}
}

/*****************************************************************************/

void
lb_cart_partition_origin(int        n_dims,
			   const int* dims,
			   const int* gsize,
			   const int* coords,
			   int*       origin)
{
	int d;

	lb_assert(n_dims > 0);

	lb_assert(dims != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(dims[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(gsize != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(gsize[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(coords != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(coords[d] >= 0 && coords[d] < dims[d]);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(origin != NULL);

	for (d = 0; d < n_dims; ++d)
		origin[d] = (gsize[d]/dims[d])*coords[d];
}

/*****************************************************************************/

void
lb_cart_partition(int        n_dims,
		    const int* gsize,
		    const int* periods,
		    int        n_processors,
		    int*       dims)
{
	int  d;
	int* n;
	int  min_cost;
	int  done;

	lb_assert(n_dims > 0);

	lb_assert(gsize != NULL);
#ifndef LB_DISABLE_ASSERTS
	for (d = 0; d < n_dims; ++d)
		lb_assert(gsize[d] > 0);
#endif /* LB_DISABLE_ASSERTS */

	lb_assert(periods != NULL);
	lb_assert(dims != NULL);

	n = lb_newa(int, n_dims);
	for (d = 1; d < n_dims; ++d)
		dims[d] = n[d] = 1;

	dims[0] = n[0] = n_processors;
	min_cost = brick_communication_cost(n_dims, n, gsize, periods);

	/*
	 *  quite stupid
	 */

	for (done = FALSE; done == FALSE;) {
		int tmp;

		for (d = 0; d < n_dims; ++d)
			if (n[d] != n_processors) {
				++n[d];
				break;
			} else {
				n[d] = 1;
				if (d == n_dims - 1)
					done = TRUE;
			}

		for (d = 0, tmp = 1; d < n_dims; ++d)
			tmp *= n[d];

		if (tmp == n_processors) {
			int cost;

			cost = brick_communication_cost
					    (n_dims, n, gsize, periods);

			if (cost < min_cost) {
				min_cost = cost;
				for (d = 0; d < n_dims; ++d)
					dims[d] = n[d];
			}
		}
	}
}
