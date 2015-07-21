#ifndef LB_CART_PARTITION_H
#define LB_CART_PARTITION_H

#include <lb-macros.h>

LB_BEGIN_DECLS

void
lb_cart_partition_size(int        n_dimensions,
			 const int* dimensions,
			 const int* global_size,
			 const int* processor_coords,
			 int*       local_size);

void
lb_cart_partition_origin(int        n_dimensions,
			   const int* dimensions,
			   const int* global_size,
			   const int* processor_coords,
			   int*       origin);

void
lb_cart_partition(int        n_dimensions,
		    const int* global_size,
		    const int* periods,
		    int        n_processors,
		    int*       dimensions);

LB_END_DECLS

#endif /* LB_CART_PARTITION_H */
