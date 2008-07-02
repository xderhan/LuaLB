#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "lb-memory.h"
#include "lb-mpi-utils.h"
#include "lb-cart-partition.h"

int
lb_mpi_comm_size(MPI_Comm comm)
{
	int size;
	LB_CALL_MPI(MPI_Comm_size(comm, &size));
	return size;
}

int
lb_mpi_comm_rank(MPI_Comm comm)
{
	int rank;
	LB_CALL_MPI(MPI_Comm_rank(comm, &rank));
	return rank;
}

static void
print_line(void)
{
	int i;

	for (i = 0; i < 64; ++i)
		putchar('=');

	putchar('\n');
}

static void
print_int_array(int n, const int* x,
		const char* start, const char* separator, const char* end)
{
	int i;

	lb_assert(n > 0);
	lb_assert(x != NULL);

	if (start != NULL)
		fputs(start, stdout);

	for (i = 0; i < n - 1; ++i)
		printf("%d%s", x[i], (separator == NULL ? ", " : separator));

	printf("%d", x[n - 1]);

	if (end != NULL)
		fputs(end, stdout);
}

static void
print_cart_info_one(int rank, const char* name, int pid,
		    int ndims,
		    const int* coords,
		    const int* size,
		    const int* origin)
{
	print_line();
	printf("<>>     process rank : %d\n"
    	       "<>>  machine : '%s', pid : %d\n",
		rank, name, pid);

	print_int_array(ndims, coords,
		        "<>>    coordinates = (", ", ", ")\n");
	print_int_array(ndims, size,
			"<>>     brick size = ", " x ", "\n");
	print_int_array(ndims, origin,
			"<>>         origin = (", ", ", ")\n");
	fflush(stdout);
}

static void
print_cart_info(MPI_Comm comm, int ndims,
	        const int* coords, const int* size, const int* origin)
{
	int rank = lb_mpi_comm_rank(comm);

	int tmp, data_size;
	char* data;

	data_size = MPI_MAX_PROCESSOR_NAME + sizeof(int) + 3*ndims*sizeof(int);
	data = lb_newa(char, data_size);

	MPI_Get_processor_name(data, &tmp);
	tmp = ndims*sizeof(int);

	if (rank) {
		int pid;

		pid = getpid();

		memcpy(data + MPI_MAX_PROCESSOR_NAME,
			&pid, sizeof(int));

		memcpy(data + MPI_MAX_PROCESSOR_NAME + sizeof(pid),
			coords, tmp);

		memcpy(data + MPI_MAX_PROCESSOR_NAME + sizeof(pid) + tmp,
			size, tmp);

		memcpy(data + MPI_MAX_PROCESSOR_NAME + sizeof(pid) + 2*tmp,
			origin, tmp);

		MPI_Send(data, data_size, MPI_CHAR, 0, 0, comm);
	} else {
	    int i;

	    print_cart_info_one(0, data, getpid(), ndims,
				coords, size, origin);

	    for (i = 1; i < lb_mpi_comm_size(comm); ++i) {
		    MPI_Status status;

		    MPI_Recv(data, data_size, MPI_CHAR, i, 0, comm, &status);

		    print_cart_info_one(i, data,
			*((int *) (data + MPI_MAX_PROCESSOR_NAME)), ndims,
			(int *) (data + MPI_MAX_PROCESSOR_NAME + sizeof(int)),
			(int *) (data + MPI_MAX_PROCESSOR_NAME + sizeof(int) + tmp),
			(int *) (data + MPI_MAX_PROCESSOR_NAME + sizeof(int) + 2*tmp));
	    }

	    print_line();
	    fflush(stdout);
	}
}

void
lb_mpi_cart_create(MPI_Comm comm, int n_dimensions,
		     const int* global_size,
		     const int* global_periods,
		     int*       partition_dims,
		     int*       processor_coords,
		     int*       local_size,
		     int*       local_origin,
		     MPI_Comm*  cart_comm)
{
	lb_assert(comm != MPI_COMM_NULL);
	lb_assert(n_dimensions > 0);
	lb_assert(global_size && global_periods);
	lb_assert(partition_dims && processor_coords);
	lb_assert(local_size && local_origin);
	lb_assert(cart_comm != NULL);

	/*
	 *  look for optimal domain partition
	 */

	lb_cart_partition(n_dimensions, global_size, global_periods,
			    lb_mpi_comm_size(comm), partition_dims);

	/*
	 *  create communicator
	 */

	LB_CALL_MPI(MPI_Cart_create(comm, n_dimensions, partition_dims,
				      (int*) global_periods, 0, cart_comm))

	/*
	 *  ask for our coordinates within processors grid
	 */

	LB_CALL_MPI(MPI_Cart_coords(*cart_comm, lb_mpi_comm_rank(*cart_comm),
				      n_dimensions, processor_coords))

	/*
	 *  setup local size/origin
	 */

	lb_cart_partition_size(n_dimensions, partition_dims, global_size,
				 processor_coords, local_size);

	lb_cart_partition_origin(n_dimensions, partition_dims, global_size,
				   processor_coords, local_origin);

	if (lb_info_enabled())
		print_cart_info(*cart_comm, n_dimensions, processor_coords,
				local_size, local_origin);
}

int
lb_mpi_neighbor_rank(MPI_Comm comm, const int* shift)
{
	int  ndims, i;
	int* dims;
	int* pers;
	int* coords;
	int  ok, rank = MPI_PROC_NULL;

	lb_assert(comm != MPI_COMM_NULL);
	lb_assert(shift != NULL);

	LB_CALL_MPI(MPI_Cartdim_get(comm, &ndims))
	lb_assert(ndims > 0);

	dims = lb_newa(int, ndims);
	pers = lb_newa(int, ndims);
	coords = lb_newa(int, ndims);

	LB_CALL_MPI(MPI_Cart_get(comm, ndims, dims, pers, coords))

	for (i = 0, ok = TRUE; i < ndims; ++i) {
		coords[i] += shift[i];
		ok = ok && (pers[i] || (coords[i] >= 0 && coords[i] < dims[i]));
	}

	if (ok)
		LB_CALL_MPI(MPI_Cart_rank(comm, coords, &rank))

	return rank;
}

void
lb_mpi_barrier(void)
{

#ifdef LB_ENABLE_MPI
	LB_CALL_MPI(MPI_Barrier(MPI_COMM_WORLD))
#endif /* LB_ENABLE_MPI */

}
