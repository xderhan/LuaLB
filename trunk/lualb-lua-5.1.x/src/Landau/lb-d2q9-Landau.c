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
#include <math.h>
#include <hdf5.h>

#include "lb-d2q9.h"
#include "lb-d2q9-utils.h"
#include "lb-memory.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_MPI
#include <mpi.h>
#include "lb-mpi-utils.h"
#endif /* LB_ENABLE_MPI */

#include "lb-d2q9-Landau.h"

struct _LBD2Q9Landau {
	LBD2Q9          d2q9;
	LBLandauParameters parameters;
	double*         pdfs[9];
	double*         rho;
	double*         u;
	double          ux_b1;
	double          ux_b3;
#ifdef LB_ENABLE_MPI
	MPI_Datatype    col_type;
	MPI_Datatype    row_type;
#endif /* LB_ENABLE_MPI */
};

/*
 *  handy macros
 */

#define NX (self->d2q9.partition_info.size[0])
#define NY (self->d2q9.partition_info.size[1])
#define RANK (self->d2q9.particle_rank)
#define OFFSET(x,y) RANK*( (x)*(NY + 2) + (y) )

/*
 *  private methods
 */

static void
lb_d2q9_Landau_average_pdfs_1(lb_d2q9_Landau* self, int o)
{
	double t1 = self->pdfs[0][o] + self->pdfs[7][o];
	double t2 = self->pdfs[1][o] + self->pdfs[5][o];
	double t3 = self->pdfs[4][o] - self->pdfs[6][o];
	
	double rho = t1 + t2 + self->pdfs[2][o]
					 + self->pdfs[3][o] + self->pdfs[4][o]
					 + self->pdfs[6][o] + self->pdfs[8][o];
	
	self->rho[o] = rho;

	self->u[o + 0] = (t1 + t3 - self->pdfs[2][o] - self->pdfs[5][o])/rho;
	self->u[o + 1] = (t2 + t3 - self->pdfs[3][o] - self->pdfs[7][o])/rho;

	o++;

	/* phi */
	t1 = self->pdfs[0][o] + self->pdfs[7][o];
	t2 = self->pdfs[1][o] + self->pdfs[5][o];
	t3 = self->pdfs[4][o] - self->pdfs[6][o];
	
	rho = t1 + t2 + self->pdfs[2][o]
					 + self->pdfs[3][o] + self->pdfs[4][o]
					 + self->pdfs[6][o] + self->pdfs[8][o];

	self->rho[o] = rho;
}

static void
lb_d2q9_Landau_average_pdfs(lb_d2q9_Landau* self)
{
	int x, y;

	lb_assert(self != NULL);

	for (x = 1; x <= NX; ++x)
		for (y = 1; y <= NY; ++y)
			lb_d2q9_Landau_average_pdfs_1(self, OFFSET(x, y) );
}

#ifdef LB_ENABLE_MPI
static void
lb_d2q9_Landau_exchange_rho(lb_d2q9_Landau* self)
{
	const int todo[4][4] = { /* Xsend, Ysend, Xrecv, Yrecv */
		{NX, 1, NX + 1, 1},
		{1, NY, 1, NY + 1},
		{1, 1, 0, 1},
		{1, 1, 1, 0}
	};

	int n, nr;

	for (n = 0, nr = 0; n < 4; ++n) {
		int rank = lb_d2q9_need_exchange(&self->d2q9, n);

		if (rank != MPI_PROC_NULL) {
			double* send = self->rho
				     + OFFSET(todo[n][0], todo[n][1]);

			double* recv = self->rho
				     + OFFSET(todo[n][2], todo[n][3]);

			MPI_Datatype type =
				(n % 2 ? self->row_type : self->col_type);

			LB_CALL_MPI(MPI_Isend(send, 1, type, rank, n,
				self->d2q9.comm, self->d2q9.requests + nr++))

			LB_CALL_MPI(MPI_Irecv(recv, 1, type, rank, (n + 2)%4,
				self->d2q9.comm, self->d2q9.requests + nr++))
		}
	}

	if (nr > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, nr);

		LB_CALL_MPI(MPI_Waitall(nr, self->d2q9.requests, statuses))
	}
}
#endif /* LB_ENABLE_MPI */

static void
lb_d2q9_Landau_wrap_rho(lb_d2q9_Landau* self)
{
	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
        if (self->d2q9.neighbors[0] ==
	    self->d2q9.partition_info.processor_rank) {
#else
        if (self->d2q9.partition_info.periods[0]) {
#endif /* LB_ENABLE_MPI */
		lb_scpy(self->rho + OFFSET(NX + 1, 1),
			    self->rho + OFFSET(1, 1),
				2, 2, NY);

		lb_scpy(self->rho + OFFSET(0, 1),
				self->rho + OFFSET(NX, 1),
				2, 2, NY);
	}

#ifdef LB_ENABLE_MPI
        if (self->d2q9.neighbors[1] ==
	    self->d2q9.partition_info.processor_rank) {
#else
        if (self->d2q9.partition_info.periods[1]) {
#endif /* LB_ENABLE_MPI */
		lb_scpy(self->rho + OFFSET(1, NY + 1),
				self->rho + OFFSET(1, 1),
				2, 2*(NY + 2), NX);

		lb_scpy(self->rho + OFFSET(1, 0),
				self->rho + OFFSET(1, NY),
				2, 2*(NY + 2), NX);
	}
}

static void
lb_d2q9_Landau_setup_rho_guards(lb_d2q9_Landau* self)
{
	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
	lb_d2q9_Landau_exchange_rho(self);
#endif /* LB_ENABLE_MPI */
	lb_d2q9_Landau_wrap_rho(self);
}

static inline void
lb_d2q9_Landau_collide_node_pdfs(lb_d2q9_Landau* self, int o,
				const double* LB_RESTRICT drx,
				const double* LB_RESTRICT dry,
				const double* LB_RESTRICT dpx,
				const double* LB_RESTRICT dpy,
				const double* LB_RESTRICT ddr,
				const double* LB_RESTRICT ddp,
				const double* LB_RESTRICT rtau,
				const double* LB_RESTRICT ptau)
{
	int n;
	double* eq[9];

	for (n = 0; n < 9; ++n)
		eq[n] = lb_new(double, 2);

	//printf("drx0=%f, dry0=%f, drx1=%f, dry1=%f, ddr0=%f, ddr1=%f\n", *drx0,*dry0,*drx1,*dry1,*ddr0,*ddr1);

	lb_d2q9_Landau_equilibrium(&self->parameters,
				  self->rho + o, self->u + o,
				  drx, dry, dpx, dpy, ddr, ddp, eq);

	/*for (n = 0; n < 9; ++n)
		printf("eq[%d][0]=%f\n", n, eq[n][1]);*/

	for (n = 0; n < 9; ++n) {
		const double tmp = self->pdfs[n][o];
		self->pdfs[n][o] = tmp + (eq[n][0] - tmp)/(*rtau);
		const double tmp1 = self->pdfs[n][o+1];
		self->pdfs[n][o+1] = tmp1 + (eq[n][1] - tmp1)/(*ptau);
	}

	for (n = 0; n < 9; ++n)
		lb_free(eq[n]);
}

static inline void
lb_d2q9_Landau_apply_bc(double* LB_RESTRICT averages,
			   double* LB_RESTRICT u,
		       double* pdfs[9], int o,
		       const double* LB_RESTRICT U,
		       int v1, int v2)
{
	static const int V[4][2] = {
		{1, 1}, {-1, 1}, {-1, -1}, {1, -1}
	};

	const double r2 = 0.5*averages[0];
	const double p2 = 0.5*averages[1];
	const double ux = u[0] - *U;
	const double uy = u[1];

	const double t1 = r2*(ux*V[v1 - 4][0] + uy*V[v1 - 4][1]);
	const double t2 = r2*(ux*V[v2 - 4][0] + uy*V[v2 - 4][1]);

	/* rho */
	pdfs[v1][o] -= t1;
	pdfs[v2][o] -= t2;

	pdfs[8][o] += t1 + t2;

	/* velocity */
	u[0] = *U;
	u[1] = 0.0;

	/* phi */
	pdfs[v1][o+1] -= t1*p2/r2;
	pdfs[v2][o+1] -= t2*p2/r2;

	pdfs[8][o+1] += (t1 + t2)*p2/r2;

}

static void
lb_d2q9_Landau_collide_pdfs(lb_d2q9_Landau* self, double rtau, double ptau)
{
	const double zero = 0.0;

	const int y0 = (lb_d2q9_is_boundary(&self->d2q9, 3) ? 2 : 1);
	const int y1 = (lb_d2q9_is_boundary(&self->d2q9, 1) ? NY - 1 : NY);

	int x, y, r;

	/*
	 *  bulk sites first
	 */

	for (x = 1; x <= NX; ++x)
		for (y = y0; y <= y1; ++y) {
			const int offset = OFFSET(x, y);
			double drx, dry, dpx, dpy, ddr, ddp;

			{
				lb_d2q9_diff_opt(&self->d2q9, x-1, y-1, self->rho+0, &drx, &dry, &ddr);
				lb_d2q9_diff_opt(&self->d2q9, x-1, y-1, self->rho+1, &dpx, &dpy, &ddp);
			}

			lb_d2q9_Landau_collide_node_pdfs(self, offset,
					    &drx, &dry, &dpx, &dpy, 
						&ddr, &ddp, &rtau, &ptau);
		}

	/*
	 *  boundary (if any) (dry = 0 corresponds to neutral wall)
	 */

	if (lb_d2q9_is_boundary(&self->d2q9, 3))
		for (x = 1; x <= NX; ++x) {
			const int offset = OFFSET(x, 1);
			double drx, dpx, ddr, ddp;

			{
				double rx_next =
					self->rho[offset + 2*(NY + 2)];
				double rx_prev =
					self->rho[offset - 2*(NY + 2)];
				double ry_next =
					self->rho[offset + 2];

				drx = 0.5*(rx_next - rx_prev);

				ddr = rx_next + rx_prev + 0.5*ry_next
				    - 2.5*self->rho[offset];

				/* Phi */
				rx_next =
					self->rho[offset+1 + 2*(NY + 2)];
				rx_prev =
					self->rho[offset+1 - 2*(NY + 2)];
				ry_next =
					self->rho[offset+1 + 2];

				dpx = 0.5*(rx_next - rx_prev);

				ddp = rx_next + rx_prev + 0.5*ry_next
				    - 2.5*self->rho[offset+1];
			}

			lb_d2q9_Landau_apply_bc(self->rho + offset,
						   self->u + offset,
					       self->pdfs, offset,
					       &self->ux_b3, 4, 5);

			lb_d2q9_Landau_collide_node_pdfs(self, offset,
			    		    &drx, &zero, &dpx, &zero, &ddr, &ddp, &rtau, &ptau);
		}

	if (lb_d2q9_is_boundary(&self->d2q9, 1))
		for (x = 1; x <= NX; ++x) {
			const int offset = OFFSET(x, NY);
			double drx, dpx, ddr, ddp;

			{
				double rx_next =
					self->rho[offset + 2*(NY + 2)];
				double rx_prev =
					self->rho[offset - 2*(NY + 2)];
				double ry_next =
					self->rho[offset + 2];

				drx = 0.5*(rx_next - rx_prev);

				ddr = rx_next + rx_prev + 0.5*ry_next
				    - 2.5*self->rho[offset];

				rx_next =
					self->rho[offset+1 + 2*(NY + 2)];
				rx_prev =
					self->rho[offset+1 - 2*(NY + 2)];
				ry_next =
					self->rho[offset+1 + 2];

				dpx = 0.5*(rx_next - rx_prev);

				ddp = rx_next + rx_prev + 0.5*ry_next
				    - 2.5*self->rho[offset+1];
			}

			lb_d2q9_Landau_apply_bc(self->rho + offset,
						   self->u + offset,
					       self->pdfs, offset,
					       &self->ux_b1, 6, 7);

			lb_d2q9_Landau_collide_node_pdfs(self, offset,
			    		    &drx, &zero, &dpx, &zero, &ddr, &ddp, &rtau, &ptau);
		}
}

/*
 *  public interface
 */

lb_d2q9_Landau*
lb_d2q9_Landau_new(int nx, int ny, int py)
{
	int n, nnodes;
	lb_d2q9_Landau* self;

	lb_assert(nx > 2 && ny > 2);

	self = lb_new(lb_d2q9_Landau, 1);

#ifdef LB_ENABLE_MPI
	lb_d2q9_init(&self->d2q9, MPI_COMM_WORLD, 2, nx, ny, TRUE, py);

	LB_CALL_MPI(MPI_Type_vector(NY, 2, 2,
				      MPI_DOUBLE, &self->col_type))
	LB_CALL_MPI(MPI_Type_commit(&self->col_type))

	LB_CALL_MPI(MPI_Type_vector(NX, 2, 2*(NY + 2),
				      MPI_DOUBLE, &self->row_type))
	LB_CALL_MPI(MPI_Type_commit(&self->row_type))
#else
	lb_d2q9_init(&self->d2q9, 2, nx, ny, TRUE, py);
#endif /* LB_ENABLE_MPI */

	self->parameters.a = -0.1;
	self->parameters.b = 0.1;
	self->parameters.K = 0.09;
	self->parameters.rtau = 0.8;
	self->parameters.ptau = 0.8;

	if (NX < 3 || NY < 3)
		lb_fatal("D2Q9Landau: size of the grid too small"
#ifdef LB_ENABLE_MPI
			   " (or too much MPI processes)"
#endif /* LB_ENABLE_MPI */
			  );

	nnodes = 2*(NX + 2)*(NY + 2);

	for (n = 0; n < 9; ++n)
		self->pdfs[n] = lb_new(double, nnodes);

	self->rho = lb_new(double, nnodes);
	self->u = lb_new(double, nnodes);

	self->ux_b1 =
	self->ux_b3 = 0.0;

#ifdef LB_ENABLE_MPI
	LB_CALL_MPI(MPI_Barrier(self->d2q9.comm))
#endif /* LB_ENABLE_MPI */

	return self;
}

void
lb_d2q9_Landau_destroy(lb_d2q9_Landau* self)
{
	int n;

	lb_assert(self != NULL);

    	lb_free(self->rho);
		lb_free(self->u);
	for (n = 0; n < 9; ++n)
		lb_free(self->pdfs[n]);

#ifdef LB_ENABLE_MPI
	LB_CALL_MPI(MPI_Type_free(&self->row_type))
	LB_CALL_MPI(MPI_Type_free(&self->col_type))
#endif /* LB_ENABLE_MPI */

	lb_d2q9_fini(&self->d2q9);

	lb_free(self);
}

void
_lb_d2q9_Landau_partition_info(const lb_d2q9_Landau* self,
			     LBPartitionInfo* pinfo)
{
	lb_assert(self != NULL);
	lb_assert(pinfo != NULL);

	*pinfo = self->d2q9.partition_info;
}

LBPartitionInfo*
lb_d2q9_Landau_partition_info(const lb_d2q9_Landau* self)
{
	LBPartitionInfo* pinfo;

	lb_assert(self != NULL);

	pinfo = lb_new(LBPartitionInfo, 1);

	_lb_d2q9_Landau_partition_info(self, pinfo);

	return pinfo;
}

void
lb_d2q9_Landau_set_parameters(lb_d2q9_Landau*             self,
			     const LBLandauParameters* parameters)
{
	lb_assert(self != NULL);
	lb_assert(parameters != NULL);
	lb_assert(parameters->a <= 0.0);
	lb_assert(parameters->b > 0.0);
	lb_assert(parameters->K > 0.0);
	lb_assert(parameters->rtau > 0.5);
	lb_assert(parameters->ptau > 0.5);

	self->parameters = *parameters;
}

void
_lb_d2q9_Landau_get_parameters(const lb_d2q9_Landau* self,
			     LBLandauParameters* parameters)
{
	lb_assert(self != NULL);
	lb_assert(parameters != NULL);

	*parameters = self->parameters;
}

LBLandauParameters*
lb_d2q9_Landau_get_parameters(const lb_d2q9_Landau* self)
{
	LBLandauParameters* params;

	lb_assert(self != NULL);

	params = lb_new(LBLandauParameters, 1);

	_lb_d2q9_Landau_get_parameters(self, params);

	return params;
}

void
lb_d2q9_Landau_set_walls_speed(lb_d2q9_Landau* self, double vx_top, double vx_down)
{
	lb_assert(self != NULL);

	if (self->d2q9.partition_info.periods[1])
		lb_warning("lb_d2q9_Landau: setting speed of the walls will"
			     " not have any effect because of periodic boundary"
			     " conditions\n");

	self->ux_b1 = vx_top;
	self->ux_b3 = vx_down;
}

void
lb_d2q9_Landau_get_walls_speed(const lb_d2q9_Landau* self,
			      double* vx_top, double* vx_down)
{
	lb_assert(self != NULL);
	lb_assert(vx_top && vx_down);

	*vx_top = self->ux_b1;
	*vx_down = self->ux_b3;
}

void
lb_d2q9_Landau_set_averages(lb_d2q9_Landau* self, int x, int y,
			   double rho, double phi, double ux, double uy)
{
	double* ptr;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(rho > 0.0);

	const int o = OFFSET(x + 1, y + 1);

	ptr = self->rho + o; 

	*ptr++ = rho;
	*ptr = phi;
	
	ptr = self->u + o; 

	*ptr++ = ux;
	*ptr   = uy;
}

void
lb_d2q9_Landau_get_averages(const lb_d2q9_Landau* self, int x, int y,
			   double* rho, double* phi, double* ux, double* uy)
{
	const double* ptr;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(rho && ux && uy);

	const int o = OFFSET(x + 1, y + 1);

	ptr = self->rho + o; 

	*rho = *ptr++;
	*phi = *ptr;

	ptr = self->u + o; 

	*ux  = *ptr++;
	*uy  = *ptr;
}

void
lb_d2q9_Landau_set_equilibrium(lb_d2q9_Landau* self)
{
	int n;
	double* ptr;

	lb_assert(self != NULL);

	//lb_d2q9_Landau_setup_rho_guards(self);
	lb_d2q9_setup_rho_guards(&self->d2q9, self->rho);

	n = RANK*(NY + 2)*(NX + 2);
	ptr = self->pdfs[0];
	for (; n > 0; --n)
		*ptr++ = 0.0;

	for (n = 1; n < 9; ++n)
		memcpy(self->pdfs[n], self->pdfs[0],
		       RANK*(NX + 2)*(NY + 2)*sizeof(double));

	lb_d2q9_Landau_collide_pdfs(self, 1.0, 1.0);	
}

void
lb_d2q9_Landau_advance(lb_d2q9_Landau* self)
{
	lb_assert(self != NULL);

	lb_d2q9_stream_pdfs(&self->d2q9, self->pdfs);

	//lb_d2q9_Landau_average_pdfs(self);
	lb_d2q9_average_pdfs(&self->d2q9, 
				self->pdfs, self->rho, self->u);

	//lb_d2q9_Landau_setup_rho_guards(self);
	lb_d2q9_setup_rho_guards(&self->d2q9, self->rho);

	lb_d2q9_Landau_collide_pdfs(self, 
				self->parameters.rtau, 
				self->parameters.ptau);
}

void
lb_d2q9_Landau_dump(const lb_d2q9_Landau* self, const char* filename)
{
	hid_t file_id;
	hid_t file_space_id;
	hid_t mem_space_id;

	hid_t rho_set_id;
	hid_t phi_set_id;

#ifdef LB_ENABLE_MPI
	hid_t xfer_plist;
#else
# define xfer_plist H5P_DEFAULT
#endif /* LB_ENABLE_MPI */

	lb_assert(self != NULL);
	lb_assert(filename != NULL);

	/*
	 *  open the file
	 */

#ifdef LB_ENABLE_MPI
	{
		hid_t acc_tpl = H5Pcreate(H5P_FILE_ACCESS);

		H5Pset_fapl_mpio(acc_tpl, self->d2q9.comm, MPI_INFO_NULL);
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC,
				    H5P_DEFAULT, acc_tpl);
		H5Pclose(acc_tpl);
	}
#else /* !LB_ENABLE_MPI */
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif /* ICHF_ENABLE_MPI */

	/*
	 *  create spaces
	 */

	{
#ifdef LB_ENABLE_MPI
		hsize_t gdims[2] = {self->d2q9.partition_info.global_size[0],
				    self->d2q9.partition_info.global_size[1]};
#else
		hsize_t gdims[2] = {NX, NY};
#endif /* LB_ENABLE_MPI */

		hsize_t ldims[2] = {NX + 2, 2*(NY + 2)};

		file_space_id = H5Screate_simple(2, gdims, NULL);
		mem_space_id = H5Screate_simple(2, ldims, NULL);
	}

	/*
	 *  select memory slab
	 */

	{
		hsize_t start[2] = {1, 2};
		hsize_t  count[2] = {NX, NY};
		hsize_t  stride[2] = {1, 2};

		H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				    start, stride, count, NULL);
	}

	/*
	 *  select slab in the file's data space
	 */

#ifdef LB_ENABLE_MPI
	{
		hsize_t start[2] = {self->d2q9.partition_info.global_origin[0],
				     self->d2q9.partition_info.global_origin[1]};
		hsize_t  count[2] = {NX, NY};

		H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET,
				    start, NULL, count, NULL);
	}
#endif /* LB_ENABLE_MPI */

	rho_set_id = H5Dcreate(file_id, "rho", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);
	phi_set_id = H5Dcreate(file_id, "phi", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);

#ifdef LB_ENABLE_MPI
	xfer_plist = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
#endif /* LB_ENABLE_MPI */

	H5Dwrite(rho_set_id, H5T_NATIVE_DOUBLE, mem_space_id,
		 file_space_id, xfer_plist, self->rho);
	H5Dwrite(phi_set_id, H5T_NATIVE_DOUBLE, mem_space_id,
		 file_space_id, xfer_plist, self->rho+1);

#ifdef LB_ENABLE_MPI
	H5Pclose(xfer_plist);
#endif /* ICHF_ENABLE_MPI */

	H5Dclose(rho_set_id);
	H5Dclose(phi_set_id);

	H5Sclose(file_space_id);
	H5Sclose(mem_space_id);
	H5Fclose(file_id);
}

void
_lb_d2q9_Landau_stats(const lb_d2q9_Landau* self, LBD2Q9LandauStats* stats)
{
	int x, y, o;
	LBD2Q9LandauStats my;

	lb_assert(self != NULL);

	my.min_density = 1.0e15;
	my.max_density = my.max_velocity =
	my.kin_energy = my.momentum_x = my.momentum_y = 0.0;

	for (x = 0, o = 2*(NY + 3); x < NX; ++x, o += 4)
		for (y = 0; y < NY; ++y) {
			const double r  = self->rho[o];
			const double ux = self->u[o++];
			const double uy = self->u[o++];

			const double u2 = ux*ux + uy*uy;
			const double u = sqrt(u2);

			my.min_density = LB_MIN(my.min_density, r);
			my.max_density = LB_MAX(my.max_density, r);
			my.max_velocity = LB_MAX(my.max_velocity, u);
			my.kin_energy += r*u2;
			my.momentum_x += r*ux;
			my.momentum_y += r*uy;
		}

#ifdef LB_ENABLE_MPI
        LB_CALL_MPI(MPI_Allreduce(&my.min_density, &stats->min_density,
		      1, MPI_DOUBLE, MPI_MIN, self->d2q9.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.max_density, &stats->max_density,
		      1, MPI_DOUBLE, MPI_MAX, self->d2q9.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.max_velocity,  &stats->max_velocity,
		      1, MPI_DOUBLE, MPI_MAX, self->d2q9.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.kin_energy, &stats->kin_energy,
		      1, MPI_DOUBLE, MPI_SUM, self->d2q9.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.momentum_x, &stats->momentum_x,
		      1, MPI_DOUBLE, MPI_SUM, self->d2q9.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.momentum_y, &stats->momentum_y,
		      1, MPI_DOUBLE, MPI_SUM, self->d2q9.comm))

	{
		const int global_N =
			self->d2q9.partition_info.global_size[0]*
			self->d2q9.partition_info.global_size[1];

		stats->kin_energy /= 2*global_N;
		stats->momentum_x /= global_N;
		stats->momentum_y /= global_N;
	}
#else /* !LB_ENABLE_MPI */
	my.kin_energy /= 2.0*NX*NY;
	my.momentum_x /= NX*NY;
	my.momentum_y /= NX*NY;
        *stats = my;
#endif /* LB_ENABLE_MPI */
}

LBD2Q9LandauStats*
lb_d2q9_Landau_stats(const lb_d2q9_Landau* self)
{
	LBD2Q9LandauStats* stat;

	lb_assert(self != NULL);

	stat = lb_new(LBD2Q9LandauStats, 1);

	_lb_d2q9_Landau_stats(self, stat);

	return stat;
}

double
lb_d2q9_Landau_mass(const lb_d2q9_Landau* self)
{
	int nx;
	double mymass = 0.0;
	const double* rho = self->rho + OFFSET(1, 1);

	lb_assert(self != NULL);

	for (nx = NX; nx > 0; --nx) {
		int ny = NY;
		for (; ny > 0; --ny) {
			mymass += *rho;
			rho += 2;
		}
		rho += 4;
	}

#ifdef LB_ENABLE_MPI
	{
		double totalmass;
		LB_CALL_MPI(MPI_Allreduce(&mymass, &totalmass, 1,
			      MPI_DOUBLE, MPI_SUM, self->d2q9.comm))
		return totalmass;
	}
#else
	return mymass;
#endif /* LB_ENABLE_MPI */
}
