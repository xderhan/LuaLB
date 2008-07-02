#ifndef LB_PARTITION_INFO_H
#define LB_PARTITION_INFO_H

#include <lb-macros.h>

LB_BEGIN_DECLS

#define LB_MAX_DIMS 3

typedef struct {
	int ndims;

	int processor_rank;
	int processors_size[LB_MAX_DIMS];
	int processor_coords[LB_MAX_DIMS];

	int size[LB_MAX_DIMS];
	int periods[LB_MAX_DIMS];
	int global_size[LB_MAX_DIMS];
	int global_origin[LB_MAX_DIMS];
} LBPartitionInfo;

LB_END_DECLS

#endif /* LB_PARTITION_INFO_H */
