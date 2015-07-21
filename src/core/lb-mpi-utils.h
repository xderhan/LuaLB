#ifndef LB_MPI_UTILS_H
#define LB_MPI_UTILS_H

#include <mpi.h>
#include <stdio.h>
#include <lb-macros.h>
#include <lb-messages.h>

#ifndef LB_DISABLE_ASSERTS
#define LB_CALL_MPI(routine)					\
	do {							\
		int retcode = routine;				\
		if (retcode != MPI_SUCCESS)			\
			lb_fatal(LB_STRLOC			\
				   " : "			\
				   LB_STRINGIFY(routine)	\
				   " failed\n");		\
	} while(0);
#else
#define LB_CALL_MPI(routine) do { (void) routine; } while(0);
#endif /* LB_DISABLE_ASSERTS */

#define LB_MPI_EXCLUSIVE(comm, rank, code)			\
	if (lb_mpi_comm_rank(comm) == rank) {			\
		code						\
	}

LB_BEGIN_DECLS

int lb_mpi_comm_size(MPI_Comm);
int lb_mpi_comm_rank(MPI_Comm);

/*
 *  handy function for use with distributed arrays
 */

void lb_mpi_cart_create(MPI_Comm, int n_dimensions,
			  const int* global_size,
			  const int* global_periods,
			  int*       partition_dims,
			  int*       processor_coords,
			  int*       local_size,
			  int*       local_origin,
			  MPI_Comm*  cart_comm);

/*
 *  yet another handy function (communicator with cartesian topology expected)
 */

int lb_mpi_neighbor_rank(MPI_Comm, const int* shift);

void lb_mpi_barrier(void);

LB_END_DECLS

#endif /* LB_MPI_UTILS_H */
