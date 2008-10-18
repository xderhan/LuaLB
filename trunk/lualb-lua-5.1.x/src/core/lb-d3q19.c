#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <string.h>

#include "lb-d3q19.h"
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
#define NZ (self->partition_info.size[2])

#define RANK (self->particle_rank)

#define OFFSET(x,y,z) (RANK*((NZ + 2)*((NY + 2)*(x) + (y)) + (z)))
#define PDFA(i,x,y,z) (pdfs[(i)] + OFFSET(x,y,z))
#define PDFA_V(v) PDFA((v)[0], (v)[1], (v)[2], (v)[3])

#ifdef LB_ENABLE_MPI
#define THIS_PROCESS (self->partition_info.processor_rank)
#endif /* LB_ENABLE_MPI */

static void
lb_d3q19_stream_inner_pdfs(const LBD3Q19* self, double* pdfs[18])
{
	int x, y;
	size_t size;

	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

	size = RANK*NZ*sizeof(double);

	for (x = 1; x <= NX; ++x)
    		for (y = 1; y <= NY; ++y) {
			const int xe = NX + 2 - x;
			const int ye = NY + 2 - y;

			memcpy(PDFA(0, xe, y, 1), PDFA(0, xe - 1, y, 1), size);
			memcpy(PDFA(1, x - 1, y, 1), PDFA(1, x, y, 1), size);
			memcpy(PDFA(2, x, ye, 1), PDFA(2, x, ye - 1, 1), size);
			memcpy(PDFA(3, x, y - 1, 1), PDFA(3, x, y, 1), size);
			memmove(PDFA(4, x, y, 2), PDFA(4, x, y, 1), size);
			memmove(PDFA(5, x, y, 0), PDFA(5, x, y, 1), size);

			memcpy(PDFA(6, xe, ye, 1), PDFA(6, xe - 1, ye - 1, 1), size);
			memcpy(PDFA(7, x - 1, ye, 1), PDFA(7, x, ye - 1, 1), size);
			memcpy(PDFA(8, xe, y - 1, 1), PDFA(8, xe - 1, y, 1), size);
			memcpy(PDFA(9, x - 1, y - 1, 1), PDFA(9, x, y, 1), size);

			memcpy(PDFA(10, xe, y, 2), PDFA(10, xe - 1, y, 1), size);
			memcpy(PDFA(11, x - 1, y, 2), PDFA(11, x, y, 1), size);
			memcpy(PDFA(12, xe, y, 0), PDFA(12, xe - 1, y, 1), size);
			memcpy(PDFA(13, x - 1, y, 0), PDFA(13, x, y, 1), size);

			memcpy(PDFA(14, x, ye, 2), PDFA(14, x, ye - 1, 1), size);
			memcpy(PDFA(15, x, y - 1, 2), PDFA(15, x, y, 1), size);
			memcpy(PDFA(16, x, ye, 0), PDFA(16, x, ye - 1, 1), size);
			memcpy(PDFA(17, x, y - 1, 0), PDFA(17, x, y, 1), size);
		}
}

static inline void
lb_d3q19_copy_X_edge(const LBD3Q19* self, double* dst, const double* src)
{
	const size_t size = RANK*sizeof(double);
	const int stride = RANK*(NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		memcpy(dst, src, size);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_copy_Y_edge(const LBD3Q19* self, double* dst, const double* src)
{
	const size_t size = RANK*sizeof(double);
	const int stride = RANK*(NZ + 2);
	int count = NY;

	for (; count > 0; --count) {
		memcpy(dst, src, size);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_copy_Z_edge(const LBD3Q19* self, double* dst, const double* src)
{
	memcpy(dst, src, RANK*NZ*sizeof(double));
}

static inline void
lb_d3q19_copy_XY_plane(const LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		lb_d3q19_copy_Y_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_copy_XZ_plane(const LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		lb_d3q19_copy_Z_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_copy_YZ_plane(const LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NZ + 2);
	int count = NY;

	for (; count > 0; --count) {
		lb_d3q19_copy_Z_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

#define COPY_YZ(di,dx,dy,dz,si,sx,sy,sz)	\
	lb_d3q19_copy_YZ_plane(self,		\
		PDFA((di),(dx),(dy),(dz)),	\
		PDFA((si),(sx),(sy),(sz)));

static void
lb_d3q19_wrap_X(const LBD3Q19* self, double* pdfs[18])
{
	COPY_YZ(0, 1, 1, 1, 0, NX + 1, 1, 1)
	COPY_YZ(6, 1, 2, 1, 6, NX + 1, 2, 1)
	COPY_YZ(8, 1, 0, 1, 8, NX + 1, 0, 1)
	COPY_YZ(10, 1, 1, 2, 10, NX + 1, 1, 2)
	COPY_YZ(12, 1, 1, 0, 12, NX + 1, 1, 0)

	COPY_YZ(1, NX, 1, 1, 1, 0, 1, 1)
	COPY_YZ(7, NX, 2, 1, 7, 0, 2, 1)
	COPY_YZ(9, NX, 0, 1, 9, 0, 0, 1)
	COPY_YZ(11, NX, 1, 2, 11, 0, 1, 2)
	COPY_YZ(13, NX, 1, 0, 13, 0, 1, 0)
}

static void
lb_d3q19_bounce_back_0(const LBD3Q19* self, double* pdfs[18])
{
	COPY_YZ(1, NX, 1, 1, 0, NX + 1, 1, 1)
	COPY_YZ(7, NX, 1, 1, 8, NX + 1, 0, 1)
	COPY_YZ(9, NX, 1, 1, 6, NX + 1, 2, 1)
	COPY_YZ(11, NX, 1, 1, 12, NX + 1, 1, 0)
	COPY_YZ(13, NX, 1, 1, 10, NX + 1, 1, 2)
}

static void
lb_d3q19_bounce_back_1(const LBD3Q19* self, double* pdfs[18])
{
	COPY_YZ(0, 1, 1, 1, 1, 0, 1, 1)
	COPY_YZ(6, 1, 1, 1, 9, 0, 0, 1)
	COPY_YZ(8, 1, 1, 1, 7, 0, 2, 1)
	COPY_YZ(10, 1, 1, 1, 13, 0, 1, 0)
	COPY_YZ(12, 1, 1, 1, 11, 0, 1, 2)
}

#define COPY_XZ(di,dx,dy,dz,si,sx,sy,sz)	\
	lb_d3q19_copy_XZ_plane(self,		\
		PDFA((di),(dx),(dy),(dz)),	\
		PDFA((si),(sx),(sy),(sz)));

static void
lb_d3q19_wrap_Y(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XZ(2, 1, 1, 1, 2, 1, NY + 1, 1)
	COPY_XZ(6, 2, 1, 1, 6, 2, NY + 1, 1)
	COPY_XZ(7, 0, 1, 1, 7, 0, NY + 1, 1)
	COPY_XZ(14, 1, 1, 2, 14, 1, NY + 1, 2)
	COPY_XZ(16, 1, 1, 0, 16, 1, NY + 1, 0)

	COPY_XZ(3, 1, NY, 1, 3, 1, 0, 1)
	COPY_XZ(8, 2, NY, 1, 8, 2, 0, 1)
	COPY_XZ(9, 0, NY, 1, 9, 0, 0, 1)
	COPY_XZ(15, 1, NY, 2, 15, 1, 0, 2)
	COPY_XZ(17, 1, NY, 0, 17, 1, 0, 0)
}

static void
lb_d3q19_bounce_back_2(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XZ(3, 1, NY, 1, 2, 1, NY + 1, 1)
	COPY_XZ(8, 1, NY, 1, 7, 0, NY + 1, 1)
	COPY_XZ(9, 1, NY, 1, 6, 2, NY + 1, 1)
	COPY_XZ(15, 1, NY, 1, 16, 1, NY + 1, 0)
	COPY_XZ(17, 1, NY, 1, 14, 1, NY + 1, 2)
}

static void
lb_d3q19_bounce_back_3(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XZ(2, 1, 1, 1, 3, 1, 0, 1)
	COPY_XZ(6, 1, 1, 1, 9, 0, 0, 1)
	COPY_XZ(7, 1, 1, 1, 8, 2, 0, 1)
	COPY_XZ(14, 1, 1, 1, 17, 1, 0, 0)
	COPY_XZ(16, 1, 1, 1, 15, 1, 0, 2)
}

#define COPY_XY(di,dx,dy,dz,si,sx,sy,sz)	\
	lb_d3q19_copy_XY_plane(self,		\
		PDFA((di),(dx),(dy),(dz)),	\
		PDFA((si),(sx),(sy),(sz)));

static void
lb_d3q19_wrap_Z(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XY(4, 1, 1, 1, 4, 1, 1, NZ + 1)
	COPY_XY(10, 2, 1, 1, 10, 2, 1, NZ + 1)
	COPY_XY(11, 0, 1, 1, 11, 0, 1, NZ + 1)
	COPY_XY(14, 1, 2, 1, 14, 1, 2, NZ + 1)
	COPY_XY(15, 1, 0, 1, 15, 1, 0, NZ + 1)

	COPY_XY(5, 1, 1, NZ, 5, 1, 1, 0)
	COPY_XY(12, 2, 1, NZ, 12, 2, 1, 0)
	COPY_XY(13, 0, 1, NZ, 13, 0, 1, 0)
	COPY_XY(16, 1, 2, NZ, 16, 1, 2, 0)
	COPY_XY(17, 1, 0, NZ, 17, 1, 0, 0)
}

static void
lb_d3q19_bounce_back_4(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XY(5, 1, 1, NZ, 4, 1, 1, NZ + 1)
	COPY_XY(12, 1, 1, NZ, 11, 0, 1, NZ + 1)
	COPY_XY(13, 1, 1, NZ, 10, 2, 1, NZ + 1)
	COPY_XY(16, 1, 1, NZ, 15, 1, 0, NZ + 1)
	COPY_XY(17, 1, 1, NZ, 14, 1, 2, NZ + 1)
}

static void
lb_d3q19_bounce_back_5(const LBD3Q19* self, double* pdfs[18])
{
	COPY_XY(4, 1, 1, 1, 5, 1, 1, 0)
	COPY_XY(10, 1, 1, 1, 13, 0, 1, 0)
	COPY_XY(11, 1, 1, 1, 12, 2, 1, 0)
	COPY_XY(14, 1, 1, 1, 17, 1, 0, 0)
	COPY_XY(15, 1, 1, 1, 16, 1, 2, 0)
}

static void
lb_d3q19_wrap_XY(const LBD3Q19* self, double* pdfs[18])
{
	lb_d3q19_copy_Z_edge(self, PDFA(6, 1, 1, 1), PDFA(6, NX + 1, NY + 1, 1));
	lb_d3q19_copy_Z_edge(self, PDFA(7, NX, 1, 1), PDFA(7, 0, NY + 1, 1));
	lb_d3q19_copy_Z_edge(self, PDFA(8, 1, NY, 1), PDFA(8, NX + 1, 0, 1));
	lb_d3q19_copy_Z_edge(self, PDFA(9, NX, NY, 1), PDFA(9, 0, 0, 1));
}

static void
lb_d3q19_wrap_XZ(const LBD3Q19* self, double* pdfs[18])
{
	lb_d3q19_copy_Y_edge(self, PDFA(10, 1, 1, 1), PDFA(10, NX + 1, 1, NZ + 1));
	lb_d3q19_copy_Y_edge(self, PDFA(11, NX, 1, 1), PDFA(11, 0, 1, NZ + 1));
	lb_d3q19_copy_Y_edge(self, PDFA(12, 1, 1, NZ), PDFA(12, NX + 1, 1, 0));
	lb_d3q19_copy_Y_edge(self, PDFA(13, NX, 1, NZ), PDFA(13, 0, 1, 0));
}

static void
lb_d3q19_wrap_YZ(const LBD3Q19* self, double* pdfs[18])
{
	lb_d3q19_copy_X_edge(self, PDFA(14, 1, 1, 1), PDFA(14, 1, NY + 1, NZ + 1));
	lb_d3q19_copy_X_edge(self, PDFA(15, 1, NY, 1), PDFA(15, 1, 0, NZ + 1));
	lb_d3q19_copy_X_edge(self, PDFA(16, 1, 1, NZ), PDFA(16, 1, NY + 1, 0));
	lb_d3q19_copy_X_edge(self, PDFA(17, 1, NY, NZ), PDFA(17, 1, 0, 0));
}

#ifdef LB_ENABLE_MPI
static void
lb_d3q19_exchange_pdfs(const LBD3Q19* self, double* pdfs[18])
{
	/*
	 *  planes first
	 */
	const struct {
		int s[5][4], r[5][4];
	} ptodo[6] = {
		{{{0, NX + 1, 1, 1}, {6, NX + 1, 2, 1}, {8, NX + 1, 0, 1},
		  {10, NX + 1, 1, 2}, {12, NX + 1, 1, 0}},
		 {{1, NX, 1, 1}, {7, NX, 2, 1}, {9, NX, 0, 1},
		  {11, NX, 1, 2}, {13, NX, 1, 0}}},

		 {{{1, 0, 1, 1}, {7, 0, 2, 1}, {9, 0, 0, 1},
		   {11, 0, 1, 2}, {13, 0, 1, 0}},
		  {{0, 1, 1, 1}, {6, 1, 2, 1}, {8, 1, 0, 1},
		   {10, 1, 1, 2}, {12, 1, 1, 0}}},

		{{{2, 1, NY + 1, 1}, {6, 2, NY + 1, 1}, {7, 0, NY + 1, 1},
		  {14, 1, NY + 1, 2}, {16, 1, NY + 1, 0}},
		 {{3, 1, NY, 1}, {8, 2, NY, 1}, {9, 0, NY, 1},
		  {15, 1, NY, 2}, {17, 1, NY, 0}}},

		{{{3, 1, 0, 1}, {8, 2, 0, 1}, {9, 0, 0, 1},
		  {15, 1, 0, 2}, {17, 1, 0, 0}},
		 {{2, 1, 1, 1},{6, 2, 1, 1}, {7, 0, 1, 1},
		  {14, 1, 1, 2}, {16, 1, 1, 0}}},

		{{{4, 1, 1, NZ + 1}, {10, 2, 1, NZ + 1}, {11, 0, 1, NZ + 1},
		  {14, 1, 2, NZ + 1}, {15, 1, 0, NZ + 1}},
		 {{5, 1, 1, NZ}, {12, 2, 1, NZ}, {13, 0, 1, NZ},
		  {16, 1, 2, NZ}, {17, 1, 0, NZ}}},

		{{{5, 1, 1, 0},{12, 2, 1, 0}, {13, 0, 1, 0},
		  {16, 1, 2, 0}, {17, 1, 0, 0}},
		 {{4, 1, 1, 1},{10, 2, 1, 1}, {11, 0, 1, 1},
		  {14, 1, 2, 1}, {15, 1, 0, 1}}}
	};

	/*
	 *  and edges
	 */

	const struct {
		int s[4], r[4];
	} etodo[12] = {
		{{6, NX + 1, NY + 1, 1}, {9, NX, NY, 1}},
		{{7, 0, NY + 1, 1}, {8, 1, NY, 1}},
		{{8, NX + 1, 0, 1}, {7, NX, 1, 1}},
		{{9, 0, 0, 1}, {6, 1, 1, 1}},

		{{10, NX + 1, 1, NZ + 1}, {13, NX, 1, NZ}},
		{{11, 0, 1, NZ + 1}, {12, 1, 1, NZ}},
		{{12, NX + 1, 1, 0}, {11, NX, 1, 1}},
		{{13, 0, 1, 0}, {10, 1, 1, 1}},

		{{14, 1, NY + 1, NZ + 1}, {17, 1, NY, NZ}},
		{{15, 1, 0, NZ + 1}, {16, 1, 1, NZ}},
		{{16, 1, NY + 1, 0}, {15, 1, NY, 1}},
		{{17, 1, 0, 0}, {14, 1, 1, 1}}
	};

	int n, r = 0;

	for (n = 0; n < 6; ++n) {
		MPI_Datatype type = MPI_DATATYPE_NULL;

		switch (n/2) {
			case 0: type = self->plane_YZ_type; break;
			case 1: type = self->plane_XZ_type; break;
			case 2: type = self->plane_XY_type; break;
			default: lb_assert(FALSE);
		}

		if (self->neighbors[n] != MPI_PROC_NULL &&
		    self->neighbors[n] != THIS_PROCESS) {
			int i;

			for (i = 0; i < 5; ++i) {
				MPI_Isend(PDFA_V(ptodo[n].s[i]), 1, type,
					  self->neighbors[n],
					  ptodo[n].s[i][0], self->comm,
					  self->requests + r++);
				MPI_Irecv(PDFA_V(ptodo[n].r[i]), 1, type,
					  self->neighbors[n],
					  ptodo[n].r[i][0], self->comm,
					  self->requests + r++);
			}
		}
	}

	for (n = 0; n < 12; ++n) {
		MPI_Datatype type = MPI_DATATYPE_NULL;

		switch (n/4) {
			case 0: type = self->edge_Z_type; break;
			case 1: type = self->edge_Y_type; break;
			case 2: type = self->edge_X_type; break;
			default: lb_assert(FALSE);
		}

		if (self->neighbors[n + 6] != MPI_PROC_NULL &&
		    self->neighbors[n + 6] != THIS_PROCESS) {
			MPI_Isend(PDFA_V(etodo[n].s), 1, type,
				  self->neighbors[n + 6],
				  etodo[n].s[0], self->comm, 
				  self->requests + r++);
			MPI_Irecv(PDFA_V(etodo[n].r), 1, type,
				  self->neighbors[n + 6],
				  etodo[n].r[0], self->comm,
				  self->requests + r++);
		}
	}

	lb_assert(r == self->n_requests);

	if (r > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, r);
		MPI_Waitall(r, self->requests, statuses);
	}
}
#endif /* LB_ENABLE_MPI */

static void
lb_d3q19_wrap_pdfs(const LBD3Q19* self, double* pdfs[18])
{
	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == THIS_PROCESS)
#else
        if (self->partition_info.periods[0])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_X(self, pdfs);
#ifndef LB_ENABLE_MPI
        else {
		lb_d3q19_bounce_back_0(self, pdfs);
		lb_d3q19_bounce_back_1(self, pdfs);
        }
#else /* LB_ENABLE_MPI */
        if (self->neighbors[0] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_0(self, pdfs);

        if (self->neighbors[1] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_1(self, pdfs);
#endif /* LB_ENABLE_MPI */


#ifdef LB_ENABLE_MPI
        if (self->neighbors[2] == THIS_PROCESS)
#else
        if (self->partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_Y(self, pdfs);
#ifndef LB_ENABLE_MPI
        else {
		lb_d3q19_bounce_back_2(self, pdfs);
		lb_d3q19_bounce_back_3(self, pdfs);
        }
#else /* LB_ENABLE_MPI */
        if (self->neighbors[2] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_2(self, pdfs);

        if (self->neighbors[3] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_3(self, pdfs);
#endif /* LB_ENABLE_MPI */

#ifdef LB_ENABLE_MPI
        if (self->neighbors[4] == THIS_PROCESS)
#else
        if (self->partition_info.periods[2])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_Z(self, pdfs);
#ifndef LB_ENABLE_MPI
        else {
		lb_d3q19_bounce_back_4(self, pdfs);
		lb_d3q19_bounce_back_5(self, pdfs);
        }
#else /* LB_ENABLE_MPI */
        if (self->neighbors[4] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_4(self, pdfs);

        if (self->neighbors[5] == MPI_PROC_NULL)
                lb_d3q19_bounce_back_5(self, pdfs);
#endif /* LB_ENABLE_MPI */


#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == THIS_PROCESS &&
	    self->neighbors[2] == THIS_PROCESS)
#else
        if (self->partition_info.periods[0] &&
	    self->partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_XY(self, pdfs);


#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == THIS_PROCESS &&
	    self->neighbors[4] == THIS_PROCESS)
#else
        if (self->partition_info.periods[0] &&
	    self->partition_info.periods[2])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_XZ(self, pdfs);


#ifdef LB_ENABLE_MPI
        if (self->neighbors[2] == THIS_PROCESS &&
	    self->neighbors[4] == THIS_PROCESS)
#else
        if (self->partition_info.periods[1] &&
	    self->partition_info.periods[2])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_YZ(self, pdfs);
}

/*
 *  public interface
 */

void
lb_d3q19_init(LBD3Q19* self,
#ifdef LB_ENABLE_MPI
	        MPI_Comm comm,
#endif /* LB_ENABLE_MPI */
	        int prank,
	        int nx, int ny, int nz,
	        int px, int py, int pz)
{
	lb_assert(self != NULL);
	lb_assert(prank > 0);
	lb_assert(nx > 0 && ny > 0 && nz > 0);

	self->partition_info.ndims = 3;

	self->partition_info.global_size[0] = nx;
	self->partition_info.global_size[1] = ny;
	self->partition_info.global_size[2] = nz;

	self->partition_info.periods[0] = px;
	self->partition_info.periods[1] = py;
	self->partition_info.periods[2] = pz;

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
	self->partition_info.processors_size[1] =
	self->partition_info.processors_size[2] = 1;

	self->partition_info.processor_coords[0] =
	self->partition_info.processor_coords[1] =
	self->partition_info.processor_coords[2] = 0;

	self->partition_info.size[0] = self->partition_info.global_size[0];
	self->partition_info.size[1] = self->partition_info.global_size[1];
	self->partition_info.size[2] = self->partition_info.global_size[2];

	self->partition_info.global_origin[0] =
	self->partition_info.global_origin[1] =
	self->partition_info.global_origin[2] = 0;
#endif /* LB_ENABLE_MPI */

	self->particle_rank = prank;

#ifdef LB_ENABLE_MPI

	/*
	 *  neighbors
	 */

	{
		const int shifts[18][3] = {
			{1, 0, 0}, {-1, 0, 0},
			{0, 1, 0}, {0, -1, 0},
			{0, 0, 1}, {0, 0, -1},

			{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
			{1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
			{0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}
		};

		int n;

		self->n_requests = 0;
		for (n = 0; n < 18; ++n) {
			self->neighbors[n] =
				lb_mpi_neighbor_rank(self->comm, shifts[n]);
			if (self->neighbors[n] != MPI_PROC_NULL &&
			    self->neighbors[n] != THIS_PROCESS)
				self->n_requests += (n < 6 ? 10 : 2);
		}

		if (self->n_requests > 0)
			self->requests =
				lb_new(MPI_Request, self->n_requests);
		else
			self->requests = NULL;
	}

	/*
	 *  register handy MPI types:
	 *  3 for edges && 3 for planes
	 */

	/* edges */
	MPI_Type_vector(NX, RANK, RANK*(NY + 2)*(NZ + 2),
			MPI_DOUBLE, &self->edge_X_type);
	MPI_Type_commit(&self->edge_X_type);

	MPI_Type_vector(NY, RANK, RANK*(NZ + 2),
			MPI_DOUBLE, &self->edge_Y_type);
	MPI_Type_commit(&self->edge_Y_type);

	MPI_Type_contiguous(NZ*RANK, MPI_DOUBLE, &self->edge_Z_type);
	MPI_Type_commit(&self->edge_Z_type);

	/* planes */
	MPI_Type_hvector(NX, 1, RANK*(NY + 2)*(NZ + 2)*sizeof(double),
			 self->edge_Y_type, &self->plane_XY_type);
	MPI_Type_commit(&self->plane_XY_type);

	MPI_Type_vector(NX, RANK*NZ, RANK*(NY + 2)*(NZ + 2),
			MPI_DOUBLE, &self->plane_XZ_type);
	MPI_Type_commit(&self->plane_XZ_type);

	MPI_Type_vector(NY, RANK*NZ, RANK*(NZ + 2),
			MPI_DOUBLE, &self->plane_YZ_type);
	MPI_Type_commit(&self->plane_YZ_type);
#endif /* LB_ENABLE_MPI */
}

void
lb_d3q19_fini(LBD3Q19* self)
{
	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
	lb_free(self->requests);

	MPI_Type_free(&self->plane_YZ_type);
	MPI_Type_free(&self->plane_XZ_type);
	MPI_Type_free(&self->plane_XY_type);

	MPI_Type_free(&self->edge_Z_type);
	MPI_Type_free(&self->edge_Y_type);
	MPI_Type_free(&self->edge_X_type);

	MPI_Comm_free(&self->comm);
#endif /* LB_ENABLE_MPI */
}

void
lb_d3q19_stream_pdfs(const LBD3Q19* self, double* pdfs[18])
{
	lb_assert(self != NULL);
	lb_assert(pdfs != NULL);

	lb_d3q19_stream_inner_pdfs(self, pdfs);
#ifdef LB_ENABLE_MPI
	lb_d3q19_exchange_pdfs(self, pdfs);
#endif /* LB_ENABLE_MPI */
	lb_d3q19_wrap_pdfs(self, pdfs);
}
