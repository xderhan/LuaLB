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

#include "lb-memory.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_MPI
#include "lb-mpi-utils.h"
#endif /* LB_ENABLE_MPI */

#include "lb-d2q9-utils.h"

/*
 *  private methods
 */

#define NX (self->partition_info.size[0])
#define NY (self->partition_info.size[1])

#define RANK (self->particle_rank)
#define OFFSET(x,y) ( RANK*( (NY + 2)*(x) + (y) ) )



static void
lb_d2q9_average_pdfs_1(LBD2Q9* self, const int o,

				const double **pdfs,

				double *rho,

				double *u)

{
	
	const double t1 = pdfs[0][o] + pdfs[7][o];
	
	const double t2 = pdfs[1][o] + pdfs[5][o];
	
	const double t3 = pdfs[4][o] - pdfs[6][o];
	
	
	const double rho1 = t1 + t2 + pdfs[2][o] + pdfs[3][o] + pdfs[4][o]
	
					 + pdfs[6][o] + pdfs[8][o];
	
	rho[o] = rho1;
	


	u[o    ] = (t1 + t3 - pdfs[2][o] - pdfs[5][o])/rho1;

	u[o + 1] = (t2 + t3 - pdfs[3][o] - pdfs[7][o])/rho1;


	/* Others component */
	
	{
		
			int ot;
		
			for (ot=1; ot<RANK; ot++) {
			
				const double s1 = pdfs[0][o + ot] + pdfs[7][o + ot];
	
			const double s2 = pdfs[1][o + ot] + pdfs[5][o + ot];
	
			const double s3 = pdfs[4][o + ot] - pdfs[6][o + ot];
				const double rho2 = s1 + s2 + pdfs[2][o + ot]
				
					 + pdfs[3][o + ot] + pdfs[4][o + ot]
									 + pdfs[6][o + ot] + pdfs[8][o + ot];

				rho[o + ot] = rho2;

	
		}
	}

}



void
lb_d2q9_average_pdfs(LBD2Q9* self,

			const double **pdfs,

			double *rho,

			double *u)

{
	
	int x, y;

	
	
	lb_assert(self != NULL);

	
	
	for (x = 1; x <= NX; ++x)
		
		for (y = 1; y <= NY; ++y)
		
			lb_d2q9_average_pdfs_1(self, OFFSET(x, y), 
											pdfs, rho, u);

}

#ifdef LB_ENABLE_MPI
static void
lb_d2q9_exchange_rho(LBD2Q9* self, double* rho)
{
	const int todo[4][4] = { /* Xsend, Ysend, Xrecv, Yrecv */
		{NX, 1, NX + 1, 1},
		{1, NY, 1, NY + 1},
		{1, 1, 0, 1},
		{1, 1, 1, 0}
	};

	int n, nr;

	for (n = 0, nr = 0; n < 4; ++n) {
		int rank = lb_d2q9_need_exchange(self, n);

		if (rank != MPI_PROC_NULL) {
			double* send = rho
				     + OFFSET(todo[n][0], todo[n][1]);

			double* recv = rho
				     + OFFSET(todo[n][2], todo[n][3]);

			MPI_Datatype type =
				(n % 2 ? self->row_type : self->col_type);

			LB_CALL_MPI(MPI_Isend(send, 1, type, rank, n,
				self->comm, self->requests + nr++))

			LB_CALL_MPI(MPI_Irecv(recv, 1, type, rank, (n + 2)%4,
				self->comm, self->requests + nr++))
		}
	}

	if (nr > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, nr);

		LB_CALL_MPI(MPI_Waitall(nr, self->requests, statuses))
	}
}
#endif /* LB_ENABLE_MPI */

static void
lb_d2q9_wrap_rho(LBD2Q9* self, double* rho)
{
	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
        if (self->neighbors[0] ==
	    self->partition_info.processor_rank) {
#else
        if (self->partition_info.periods[0]) {
#endif /* LB_ENABLE_MPI */
		lb_scpy(rho + OFFSET(NX + 1, 1),
			rho + OFFSET(1, 1),
				RANK, RANK, NY);

		lb_scpy(rho + OFFSET(0, 1),
			rho + OFFSET(NX, 1),
				RANK, RANK, NY);
	}

#ifdef LB_ENABLE_MPI
        if (self->neighbors[1] ==
	    self->partition_info.processor_rank) {
#else
        if (self->partition_info.periods[1]) {
#endif /* LB_ENABLE_MPI */
		lb_scpy(rho + OFFSET(1, NY + 1),
			rho + OFFSET(1, 1),
			RANK, RANK*(NY + 2), NX);

		lb_scpy(rho + OFFSET(1, 0),
			rho + OFFSET(1, NY),
			RANK, RANK*(NY + 2), NX);
	}
}

void
lb_d2q9_setup_rho_guards(LBD2Q9* self, double* rho)
{
	lb_assert(self != NULL);
	lb_assert(rho != NULL);

#ifdef LB_ENABLE_MPI
	lb_d2q9_exchange_rho(self, rho);
#endif /* LB_ENABLE_MPI */
	lb_d2q9_wrap_rho(self, rho);
}

void lb_d2q9_diff_opt(LBD2Q9* self, int x, int y, 
			const double *rho, double *drx, double *dry, double *ddr)
{
	lb_assert(self != NULL);
	lb_assert(rho != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	
	const int offset = OFFSET(x + 1, y + 1);

	/* extract for rho */
	double rx_next = rho[offset + RANK*(NY + 2)];
	double rx_prev = rho[offset - RANK*(NY + 2)];
		
	double ry_next = rho[offset + RANK];
	double ry_prev = rho[offset - RANK];

	*drx = 0.5*(rx_next - rx_prev);
	*dry = 0.5*(ry_next - ry_prev);

	*ddr = rx_next + rx_prev + ry_next + ry_prev - 4.*rho[offset];
}

void lb_d2q9_diff_opt_bb(LBD2Q9* self, int x, int y, 
			const double *rho, double *drx, double *dry, double *ddr)
{
	lb_assert(self != NULL);
	lb_assert(rho != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	
	const int offset = OFFSET(x, y);

	/* extract for rho */
	double rx_next = rho[offset + RANK*(NY + 2)];
	double rx_prev = rho[offset - RANK*(NY + 2)];
		
	double ry_next = rho[offset + RANK];
	double ry_prev = rho[offset - RANK];

	*drx = 0.5*(rx_next - rx_prev);
	*dry = 0.5*(ry_next - ry_prev);

	*ddr = rx_next + rx_prev + ry_next + ry_prev - 4.*rho[offset];
}
