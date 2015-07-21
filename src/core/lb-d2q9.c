#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <string.h>

#include "lb-d2q9.h"
#include "lb-memory.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_MPI
#include "lb-mpi-utils.h"
#endif /* LB_ENABLE_MPI */

/*
 *  private methods
 */

#define NX (self->partition_info.size[0])
#define NY (self->partition_info.size[1])

#define RANK (self->particle_rank)
#define OFFSET(x,y) (RANK*((NY + 2)*(x) + (y)))

static void
lb_d2q9_stream_inner_pdfs(const LBD2Q9* self, double* pdfs[8])
{
	int colsize, x;

	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

	colsize = RANK*NY*sizeof(double);

	for (x = 1; x <= NX; ++x) {
		memcpy(pdfs[0] + OFFSET(NX - x + 2, 1),
		       pdfs[0] + OFFSET(NX - x + 1, 1),
		       colsize);

		memmove(pdfs[1] + OFFSET(x, 2),
				pdfs[1] + OFFSET(x, 1),
		        colsize);

		memcpy(pdfs[2] + OFFSET(x - 1, 1),
		       pdfs[2] + OFFSET(x, 1),
		       colsize);

		memmove(pdfs[3] + OFFSET(x, 0),
				pdfs[3] + OFFSET(x, 1),
		        colsize);

		memcpy(pdfs[4] + OFFSET(NX - x + 2, 2),
		       pdfs[4] + OFFSET(NX - x + 1, 1),
		       colsize);

		memcpy(pdfs[5] + OFFSET(x - 1, 2),
		       pdfs[5] + OFFSET(x, 1),
		       colsize);

		memcpy(pdfs[6] + OFFSET(x - 1, 0),
		       pdfs[6] + OFFSET(x, 1),
		       colsize);

		memcpy(pdfs[7] + OFFSET(NX - x + 2, 0),
		       pdfs[7] + OFFSET(NX - x + 1, 1),
		       colsize);
	}
}

static inline double*
lb_d2q9_address(const LBD2Q9* self, double* pdfs[8], const int info[3])
{
	return pdfs[info[2]] + OFFSET(info[0], info[1]);
}

static inline void
lb_d2q9_colscpy(const LBD2Q9* self, double* pdfs[8],
		  const int info[][6], int ncols)
{
	const int colsize = RANK*NY*sizeof(double);

	for (; ncols > 0; --ncols, ++info)
		memcpy(lb_d2q9_address(self, pdfs, *info),
		       lb_d2q9_address(self, pdfs, *info + 3),
		       colsize);
}

static inline void
lb_d2q9_rowscpy(const LBD2Q9* self, double* pdfs[8],
		  const int info[][6], int nrows)
{
	for (; nrows > 0; --nrows, ++info)
		lb_scpy(lb_d2q9_address(self, pdfs, *info),
			  lb_d2q9_address(self, pdfs, *info + 3),
			  RANK, RANK*(NY + 2), NX);
}

static void
lb_d2q9_bounce_back_0(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[3][6] = {
		{NX, 1, 2, NX + 1, 1, 0},
		{NX, 1, 5, NX + 1, 0, 7},
		{NX, 1, 6, NX + 1, 2, 4}
	};

	lb_d2q9_colscpy(self, pdfs, todo, 3);
}

static void
lb_d2q9_bounce_back_1(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[3][6] = {
		{1, NY, 3, 1, NY + 1, 1},
		{1, NY, 6, 2, NY + 1, 4},
		{1, NY, 7, 0, NY + 1, 5}
	};

	lb_d2q9_rowscpy(self, pdfs, todo, 3);
}

static void
lb_d2q9_bounce_back_2(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[3][6] = {
		{1, 1, 0, 0, 1, 2},
		{1, 1, 4, 0, 0, 6},
		{1, 1, 7, 0, 2, 5}
	};

	lb_d2q9_colscpy(self, pdfs, todo, 3);
}

static void
lb_d2q9_bounce_back_3(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[3][6] = {
		{1, 1, 1, 1, 0, 3},
		{1, 1, 4, 0, 0, 6},
		{1, 1, 5, 2, 0, 7}
	};

	lb_d2q9_rowscpy(self, pdfs, todo, 3);
}

static void
lb_d2q9_wrap_0_2(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[6][6] = {
		{NX, 1, 2, 0, 1, 2},
		{NX, 2, 5, 0, 2, 5},
		{NX, 0, 6, 0, 0, 6},
		{1, 1, 0, NX + 1, 1, 0},
		{1, 2, 4, NX + 1, 2, 4},
		{1, 0, 7, NX + 1, 0, 7}
	};

	lb_d2q9_colscpy(self, pdfs, todo, 6);
}

static void
lb_d2q9_wrap_1_3(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[6][6] = {
		{0, NY, 6, 0, 0, 6},
		{1, NY, 3, 1, 0, 3},
		{2, NY, 7, 2, 0, 7},
		{0, 1, 5, 0, NY + 1, 5},
		{1, 1, 1, 1, NY + 1, 1},
		{2, 1, 4, 2, NY + 1, 4}
	};

	lb_d2q9_rowscpy(self, pdfs, todo, 6);
}

static void
lb_d2q9_wrap_corners(const LBD2Q9* self, double* pdfs[8])
{
	const int todo[4][6] = {
		{NX, NY, 6, 0, 0, 6},
		{NX, 1, 5, 0, NY + 1, 5},
		{1, NY, 7, NX + 1, 0, 7},
		{1, 1, 4, NX + 1, NY + 1, 4}
	};

	int i;

	for (i = 0; i < 4; ++i)
		memcpy(pdfs[todo[i][2]] + OFFSET(todo[i][0], todo[i][1]),
		       pdfs[todo[i][5]] + OFFSET(todo[i][3], todo[i][4]),
		       RANK*sizeof(double));
}

static void
lb_d2q9_wrap_pdfs(const LBD2Q9*self, double* pdfs[8])
{
	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == self->partition_info.processor_rank)
#else
        if (self->partition_info.periods[0])
#endif /* LB_ENABLE_MPI */
		lb_d2q9_wrap_0_2(self, pdfs);
#ifndef LB_ENABLE_MPI
        else {
		lb_d2q9_bounce_back_0(self, pdfs);
		lb_d2q9_bounce_back_2(self, pdfs);
        }
#else /* LB_ENABLE_MPI */
        if (self->neighbors[0] == MPI_PROC_NULL)
                lb_d2q9_bounce_back_0(self, pdfs);

        if (self->neighbors[2] == MPI_PROC_NULL)
                lb_d2q9_bounce_back_2(self, pdfs);
#endif /* LB_ENABLE_MPI */

#ifdef LB_ENABLE_MPI
        if (self->neighbors[1] == self->partition_info.processor_rank)
#else
        if (self->partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d2q9_wrap_1_3(self, pdfs);
#ifndef LB_ENABLE_MPI
        else {
		lb_d2q9_bounce_back_1(self, pdfs);
		lb_d2q9_bounce_back_3(self, pdfs);
        }
#else /* LB_ENABLE_MPI */
        if (self->neighbors[1] == MPI_PROC_NULL)
                lb_d2q9_bounce_back_1(self, pdfs);

        if (self->neighbors[3] == MPI_PROC_NULL)
                lb_d2q9_bounce_back_3(self, pdfs);
#endif /* LB_ENABLE_MPI */

#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == self->partition_info.processor_rank &&
            self->neighbors[1] == self->partition_info.processor_rank)
#else
        if (self->partition_info.periods[0] &&
	    self->partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d2q9_wrap_corners(self, pdfs);
}

#ifdef LB_ENABLE_MPI
static void
lb_d2q9_exchange_pdfs(const LBD2Q9* self, double* pdfs[8])
{
	const int exch1[4][6][3] = {
		{{NX + 1, 1, 0}, {NX, 1, 2},
		 {NX + 1, 2, 4}, {NX, 2, 5},
		 {NX + 1, 0, 7}, {NX, 0, 6}},

		{{1, NY + 1, 1}, {1, NY, 3},
		 {2, NY + 1, 4}, {2, NY, 7},
		 {0, NY + 1, 5}, {0, NY, 6}},

		{{0, 1, 2}, {1, 1, 0}, {0, 2, 5},
		 {1, 2, 4}, {0, 0, 6}, {1, 0, 7}},

		{{1, 0, 3}, {1, 1, 1}, {2, 0, 7},
		 {2, 1, 4}, {0, 0, 6}, {0, 1, 5}}
	};

	const int exch2[4][2][3] = {
		{{NX + 1, NY + 1, 4}, {NX, NY, 6}},
		{{0, NY + 1, 5}, {1, NY, 7}},
		{{0, 0, 6}, {1, 1, 4}},
		{{NX + 1, 0, 7}, {NX, 1, 5}}
	};

	int n, i, r;

	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

	for (n = 0, r = 0; n < 4; ++n)
		if (lb_d2q9_need_exchange(self, n) != MPI_PROC_NULL) {
			MPI_Datatype type =
				(n % 2 ? self->row_type : self->col_type);

			for (i = 0; i < 6; i += 2) {
				double* send = pdfs[exch1[n][i][2]]
				 + OFFSET(exch1[n][i][0], exch1[n][i][1]);

				double* recv = pdfs[exch1[n][i + 1][2]]
				 + OFFSET(exch1[n][i + 1][0], exch1[n][i + 1][1]);

				LB_CALL_MPI(MPI_Isend(send, 1, type,
				    	      self->neighbors[n],
					      exch1[n][i][2], self->comm,
				    	      self->requests + r++))

				LB_CALL_MPI(MPI_Irecv(recv, 1, type,
				    	      self->neighbors[n],
				    	      exch1[n][i + 1][2], self->comm,
				    	      self->requests + r++))
			}
		}

	for (n = 0; n < 4; ++n)
		if (lb_d2q9_need_exchange(self, n + 4) != MPI_PROC_NULL) {
			double* send = pdfs[exch2[n][0][2]]
				     + OFFSET(exch2[n][0][0], exch2[n][0][1]);

			double* recv = pdfs[exch2[n][1][2]]
				     + OFFSET(exch2[n][1][0], exch2[n][1][1]);

			LB_CALL_MPI(MPI_Isend(send, self->particle_rank,
				      MPI_DOUBLE, self->neighbors[n + 4],
				      exch2[n][0][2], self->comm,
				      self->requests + r++))

			LB_CALL_MPI(MPI_Irecv(recv, self->particle_rank,
				      MPI_DOUBLE, self->neighbors[n + 4],
				      exch2[n][1][2], self->comm,
				      self->requests + r++))
		}

	if (self->n_requests > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, r);

		lb_assert(r == self->n_requests);

		LB_CALL_MPI(MPI_Waitall(r, self->requests, statuses))
	}
}
#endif /* LB_ENABLE_MPI */

/*
 *  public interface
 */

void
lb_d2q9_init(LBD2Q9* self,
#ifdef LB_ENABLE_MPI
	       MPI_Comm comm,
#endif /* LB_ENABLE_MPI */
	       int prank,
	       int nx, int ny,
	       int px, int py)
{
	lb_assert(self != NULL);
	lb_assert(prank > 0);
	lb_assert(nx > 0 && ny > 0);

	self->partition_info.ndims = 2;

	self->partition_info.global_size[0] = nx;
	self->partition_info.global_size[1] = ny;
	self->partition_info.periods[0] = px;
	self->partition_info.periods[1] = py;

#ifdef LB_ENABLE_MPI
	lb_mpi_cart_create(comm,
			     self->partition_info.ndims,
			     self->partition_info.global_size,
			     self->partition_info.periods,
			     self->partition_info.processors_size,
			     self->partition_info.processor_coords,
			     self->partition_info.size,
			     self->partition_info.global_origin,
			     &self->comm);

	self->partition_info.processor_rank =
		lb_mpi_comm_rank(self->comm);

#else /* !LB_ENABLE_MPI */
	self->partition_info.processor_rank = 0;

	self->partition_info.processors_size[0] =
	self->partition_info.processors_size[1] = 1;

	self->partition_info.processor_coords[0] =
	self->partition_info.processor_coords[1] = 0;

	self->partition_info.size[0] = self->partition_info.global_size[0];
	self->partition_info.size[1] = self->partition_info.global_size[1];

	self->partition_info.global_origin[0] =
	self->partition_info.global_origin[1] = 0;
#endif /* LB_ENABLE_MPI */

	self->particle_rank = prank;

#ifdef LB_ENABLE_MPI
	/*
	 *  neighbors
	 */

	{
		const int shifts[8][2] = {
			{1, 0}, { 0, 1}, {-1,  0}, {0, -1},
			{1, 1}, {-1, 1}, {-1, -1}, {1, -1},
		};

		int n;

		for (n = 0; n < 8; ++n)
			self->neighbors[n] =
				lb_mpi_neighbor_rank(self->comm, shifts[n]);

		self->n_requests = 0;
		for (n = 0; n < 4; ++n)
			if (lb_d2q9_need_exchange(self, n) != MPI_PROC_NULL)
				self->n_requests += 6;

		for (n = 4; n < 8; ++n)
			if (lb_d2q9_need_exchange(self, n) != MPI_PROC_NULL)
				self->n_requests += 2;

		if (self->n_requests > 0)
			self->requests =
				lb_new(MPI_Request, self->n_requests);
		else
			self->requests = NULL;
	}

	/*
	 *  register handy MPI types
	 */

	LB_CALL_MPI(MPI_Type_contiguous(NY*RANK, MPI_DOUBLE, &self->col_type))
	LB_CALL_MPI(MPI_Type_commit(&self->col_type))

	LB_CALL_MPI(MPI_Type_vector(NX, RANK, RANK*(NY + 2),
				      MPI_DOUBLE, &self->row_type))
	LB_CALL_MPI(MPI_Type_commit(&self->row_type))
#endif /* LB_ENABLE_MPI */
}

void
lb_d2q9_fini(LBD2Q9* self)
{
	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
	lb_free(self->requests);

	LB_CALL_MPI(MPI_Type_free(&self->row_type))
	LB_CALL_MPI(MPI_Type_free(&self->col_type))
	LB_CALL_MPI(MPI_Comm_free(&self->comm))
#endif /* LB_ENABLE_MPI */
}

void
lb_d2q9_stream_pdfs(const LBD2Q9* self, double* pdfs[8])
{
	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

	lb_d2q9_stream_inner_pdfs(self, pdfs);
#ifdef LB_ENABLE_MPI
	lb_d2q9_exchange_pdfs(self, pdfs);
#endif /* LB_ENABLE_MPI */
	lb_d2q9_wrap_pdfs(self, pdfs);
}

#ifdef LB_ENABLE_MPI
int
lb_d2q9_need_exchange(const LBD2Q9* self, int n)
{
	lb_assert(self != NULL);
	lb_assert(n >= 0 && n <= 7);

	if (self->neighbors[n] != MPI_PROC_NULL &&
	    self->neighbors[n] != self->partition_info.processor_rank)
		return self->neighbors[n];
	else
		return MPI_PROC_NULL;
}
#endif /* LB_ENABLE_MPI */

int
lb_d2q9_is_boundary(const LBD2Q9* self, int n)
{
	lb_assert(self != NULL);
	lb_assert(n >= 0 && n <= 7);

#ifdef LB_ENABLE_MPI
	return self->neighbors[n] == MPI_PROC_NULL;
#else
	return n < 4 ? !self->partition_info.periods[n % 2]
		     : (!self->partition_info.periods[0] ||
		        !self->partition_info.periods[1]);
#endif /* LB_ENABLE_MPI */
}
