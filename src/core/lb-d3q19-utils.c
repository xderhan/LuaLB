/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <string.h>

#include "lb-d3q19-utils.h"
#include "lb-messages.h"
#include "lb-memory.h"

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
#define OFFSET_V(a) (OFFSET((a)[0], (a)[1], (a)[2]))

#ifdef LB_ENABLE_MPI
#define THIS_PROCESS (self->partition_info.processor_rank)
#endif /* LB_ENABLE_MPI */



static inline void
lb_d3q19_average_pdfs_1(LBD3Q19* self, const int o,

				const double **pdfs,

				double *rho,

				double *u)

{
	
	const double t1 = pdfs[6][o]  + pdfs[8][o]
	
				+ pdfs[10][o] + pdfs[12][o];
	const double t2 = pdfs[7][o]  + pdfs[9][o]
	
		        + pdfs[11][o] + pdfs[13][o];
	
	const double t3 = pdfs[6][o]  + pdfs[7][o]
		
			+ pdfs[14][o] + pdfs[16][o];
	
	const double t4 = pdfs[8][o]  + pdfs[9][o]
		
			+ pdfs[15][o] + pdfs[17][o];
	
	const double t5 = pdfs[10][o] + pdfs[11][o]
		
			+ pdfs[14][o] + pdfs[15][o];
	
	const double t6 = pdfs[12][o] + pdfs[13][o]
		
			+ pdfs[16][o] + pdfs[17][o];
	
	const double rho1 = pdfs[0][o] + pdfs[1][o]
		
			  + pdfs[2][o] + pdfs[3][o]
		
			  + pdfs[4][o] + pdfs[5][o]
		
			  + 0.5*(t1 + t2 + t3 + t4 + t5 + t6)
		
			  + pdfs[18][o];
	
	rho[o] = rho1;
	
	u[o    ] = (pdfs[0][o] + t1 - pdfs[1][o] - t2)/rho1;
	u[o + 1] = (pdfs[2][o] + t3 - pdfs[3][o] - t4)/rho1;
	u[o + 2] = (pdfs[4][o] + t5 - pdfs[5][o] - t6)/rho1;
	
	/* For Others Component */
	
	{
		
		int ot;
		
		for (ot=1; ot<RANK; ot++) {
			
			const double s1 = pdfs[6][o + ot]  + pdfs[8][o + ot]
			
					+ pdfs[10][o + ot] + pdfs[12][o + ot];

			

			const double s2 = pdfs[7][o + ot]  + pdfs[9][o + ot]
			
		        		+ pdfs[11][o + ot] + pdfs[13][o + ot];
			
			

			const double s3 = pdfs[6][o + ot]  + pdfs[7][o + ot]
			
					+ pdfs[14][o + ot] + pdfs[16][o + ot];
			
			

			const double s4 = pdfs[8][o + ot]  + pdfs[9][o + ot]
			
					+ pdfs[15][o + ot] + pdfs[17][o + ot];
			
			

			const double s5 = pdfs[10][o + ot] + pdfs[11][o + ot]
			
					+ pdfs[14][o + ot] + pdfs[15][o + ot];
			
			

			const double s6 = pdfs[12][o + ot] + pdfs[13][o + ot]
			
					+ pdfs[16][o + ot] + pdfs[17][o + ot];
			
			

			const double rho2 = pdfs[0][o + ot] + pdfs[1][o + ot]
			
			 		+ pdfs[2][o + ot] + pdfs[3][o + ot]
				
			 		+ pdfs[4][o + ot] + pdfs[5][o + ot]
				
			 		+ 0.5*(s1 + s2 + s3 + s4 + s5 + s6)
				
			 		+ pdfs[18][o + ot];
			
			

			rho[o + ot] = rho2;
			
			
		}
		}
}

void

lb_d3q19_average_pdfs(LBD3Q19* self,

			const double **pdfs,

			double *rho,

			double *u)
{
	
	int x, y, z;
	
	

	lb_assert(self != NULL);
	
	

	for (x = 1; x <= NX; ++x)
	
		for (y = 1; y <= NY; ++y)
		
			for (z = 1; z <= NZ; ++z)
		
				lb_d3q19_average_pdfs_1
		
					(self, OFFSET(x, y, z), pdfs, rho, u);
}


#ifdef LB_ENABLE_MPI
static void
lb_d3q19_exchange_rho(LBD3Q19* self, double *rho)
{
	static const int opposite[6] = {1, 0, 3, 2, 5, 4};

	const struct {
		int s[3], r[3];
	} todo[6] = {
		{{NX, 1, 1}, {NX + 1, 1, 1}},
		{{1, 1, 1}, {0, 1, 1}},
		{{1, NY, 1}, {1, NY + 1, 1}},
		{{1, 1, 1}, {1, 0, 1}},
		{{1, 1, NZ}, {1, 1, NZ + 1}},
		{{1, 1, 1}, {1, 1, 0}}
	};

	int n, r = 0;

	for (n = 0; n < 6; ++n)
		if (self->neighbors[n] != MPI_PROC_NULL &&
		    self->neighbors[n] != THIS_PROCESS) {
			MPI_Datatype type = MPI_DATATYPE_NULL;

			switch (n/2) {
				case 0 : type = self->plane_YZ_type; break;
				case 1 : type = self->plane_XZ_type; break;
				case 2 : type = self->plane_XY_type; break;
				default : lb_assert(FALSE);
			}

			MPI_Isend(rho + OFFSET_V(todo[n].s), 1, type,
				  self->neighbors[n], n,
				  self->comm, self->requests + r++);

			MPI_Irecv(rho + OFFSET_V(todo[n].r), 1, type,
				  self->neighbors[n], opposite[n],
				  self->comm, self->requests + r++);
		}

	lb_assert(r <= self->n_requests);

	if (r > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, r);
		MPI_Waitall(r, self->requests, statuses);
	}
}
#endif /* LB_ENABLE_MPI */

static inline void
lb_d3q19_copy_X_edge(LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NY + 2)*(NZ + 2);
	int count = NX;
	
	lb_scpy(dst, src, RANK, stride, count);	

	/*for (; count > 0; --count) {
		*dst = *src;
		dst += stride;
		src += stride;
	}*/
}

static inline void
lb_d3q19_copy_Y_edge(LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NZ + 2);
	int count = NY;

	lb_scpy(dst, src, RANK, stride, count);	

	/*for (; count > 0; --count) {
		*dst = *src;
		dst += stride;
		src += stride;
	}*/
}

static inline void
lb_d3q19_copy_Z_edge(LBD3Q19* self, double* dst, const double* src)
{
	memcpy(dst, src, NZ*RANK*sizeof(double));
}

static inline void
lb_d3q19_copy_XY_plane(LBD3Q19* self, double* dst, const double* src)
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
lb_d3q19_copy_XZ_plane(LBD3Q19* self, double* dst, const double* src)
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
lb_d3q19_copy_YZ_plane(LBD3Q19* self, double* dst, const double* src)
{
	const int stride = RANK*(NZ + 2);
	int count = NY;

	for (; count > 0; --count) {
		lb_d3q19_copy_Z_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static void
lb_d3q19_wrap_rho_X(LBD3Q19* self, double* rho)
{
	lb_d3q19_copy_YZ_plane(self, rho + OFFSET(0, 1, 1),
				rho + OFFSET(NX, 1, 1));
	lb_d3q19_copy_YZ_plane(self, rho + OFFSET(NX + 1, 1, 1),
				rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_wrap_rho_Y(LBD3Q19* self, double* rho)
{
	lb_d3q19_copy_XZ_plane(self, rho + OFFSET(1, 0, 1),
				rho + OFFSET(1, NY, 1));
	lb_d3q19_copy_XZ_plane(self, rho + OFFSET(1, NY + 1, 1),
				rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_wrap_rho_Z(LBD3Q19* self, double* rho)
{
	lb_d3q19_copy_XY_plane(self, rho + OFFSET(1, 1, 0),
				rho + OFFSET(1, 1, NZ));
	lb_d3q19_copy_XY_plane(self, rho + OFFSET(1, 1, NZ + 1),
				rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_wrap_rho(LBD3Q19* self, double* rho)
{
#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] == THIS_PROCESS)
#else
        if (self->partition_info.periods[0])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_rho_X(self, rho);

#ifdef LB_ENABLE_MPI
        if (self->neighbors[2] == THIS_PROCESS)
#else
        if (self->partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_rho_Y(self, rho);

#ifdef LB_ENABLE_MPI
        if (self->neighbors[4] == THIS_PROCESS)
#else
        if (self->partition_info.periods[2])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_wrap_rho_Z(self, rho);
}

void
lb_d3q19_setup_rho_guards(LBD3Q19* self, double* rho)
{
	lb_assert(self != NULL);
	lb_assert(rho != NULL);

#ifdef LB_ENABLE_MPI
	lb_d3q19_exchange_rho(self, rho);
#endif /* LB_ENABLE_MPI */
	lb_d3q19_wrap_rho(self, rho);
}
