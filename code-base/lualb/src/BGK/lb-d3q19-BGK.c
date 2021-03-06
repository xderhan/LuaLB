#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <hdf5.h>

#include "lb-d3q19.h"
#include "lb-memory.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_MPI
#include <mpi.h>
#include "lb-mpi-utils.h"
#endif /* LB_ENABLE_MPI */

#include "lb-d3q19-BGK.h"

struct _LBD3Q19BGK {
	LBD3Q19	d3q19;
	LBGKParameters parameters;
	double*	pdfs[19];
	double*	rho;
	double*	u;
	double	uxuy_b4[2];
	double	uxuy_b5[2];
};

/*
 *  handy macros
 */

#define NX (self->d3q19.partition_info.size[0])
#define NY (self->d3q19.partition_info.size[1])
#define NZ (self->d3q19.partition_info.size[2])

#define OFFSET(x,y,z) ((NZ + 2)*((x)*(NY + 2) + (y)) + (z))
#define OFFSET_V(a) (OFFSET((a)[0], (a)[1], (a)[2]))

#ifdef LB_ENABLE_MPI
#define THIS_PROCESS (self->d3q19.partition_info.processor_rank)
#endif /* LB_ENABLE_MPI */

/*
 *  private methods
 */

static inline void
lb_d3q19_BGK_average_pdfs_1(LBD3Q19BGK* self, int o)
{
	const double t1 = self->pdfs[6][o] + self->pdfs[8][o]
			+ self->pdfs[10][o] + self->pdfs[12][o];

	const double t2 = self->pdfs[7][o] + self->pdfs[9][o]
		        + self->pdfs[11][o] + self->pdfs[13][o];

	const double t3 = self->pdfs[6][o] + self->pdfs[7][o]
			+ self->pdfs[14][o] + self->pdfs[16][o];

	const double t4 = self->pdfs[8][o] + self->pdfs[9][o]
			+ self->pdfs[15][o] + self->pdfs[17][o];

	const double t5 = self->pdfs[10][o] + self->pdfs[11][o]
			+ self->pdfs[14][o] + self->pdfs[15][o];

	const double t6 = self->pdfs[12][o] + self->pdfs[13][o]
			+ self->pdfs[16][o] + self->pdfs[17][o];

	const double rho = self->pdfs[0][o] + self->pdfs[1][o]
			 + self->pdfs[2][o] + self->pdfs[3][o]
			 + self->pdfs[4][o] + self->pdfs[5][o]
			 + 0.5*(t1 + t2 + t3 + t4 + t5 + t6)
			 + self->pdfs[18][o];

	const int uo = 3*o;

	self->rho[o] = rho;

	self->u[uo + 0] = (self->pdfs[0][o] + t1 - self->pdfs[1][o] - t2)/rho;
	self->u[uo + 1] = (self->pdfs[2][o] + t3 - self->pdfs[3][o] - t4)/rho;
	self->u[uo + 2] = (self->pdfs[4][o] + t5 - self->pdfs[5][o] - t6)/rho;
}

static void
lb_d3q19_BGK_average_pdfs(LBD3Q19BGK* self)
{
	int x, y, z;

	lb_assert(self != NULL);

	for (x = 1; x <= NX; ++x)
		for (y = 1; y <= NY; ++y)
			for (z = 1; z <= NZ; ++z)
				lb_d3q19_BGK_average_pdfs_1
					(self, OFFSET(x, y, z));
}

static inline void
lb_d3q19_BGK_apply_bc(LBD3Q19BGK* self, int o,
		        const double* LB_RESTRICT U,
		        int v1, int v2, int v3, int v4)
{
	static const int V[12][3] = {
		{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
		{1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
		{0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}
	};

	const double r4 = 0.25*self->rho[o];
	const int o3 = 3*o;
	const double ux = self->u[o3 + 0] - U[0];
	const double uy = self->u[o3 + 1] - U[1];
	const double uz = self->u[o3 + 2];

	const double t1 = r4*(ux*V[v1 - 6][0] + uy*V[v1 - 6][1] + uz*V[v1 - 6][2]);
	const double t2 = r4*(ux*V[v2 - 6][0] + uy*V[v2 - 6][1] + uz*V[v2 - 6][2]);
	const double t3 = r4*(ux*V[v3 - 6][0] + uy*V[v3 - 6][1] + uz*V[v3 - 6][2]);
	const double t4 = r4*(ux*V[v4 - 6][0] + uy*V[v4 - 6][1] + uz*V[v4 - 6][2]);

	self->pdfs[v1][o] -= t1;
	self->pdfs[v2][o] -= t2;
	self->pdfs[v3][o] -= t3;
	self->pdfs[v4][o] -= t4;

	self->pdfs[18][o] += t1 + t2 + t3 + t4;

	self->u[o3 + 0] = U[0];
	self->u[o3 + 1] = U[1];
	self->u[o3 + 2] = 0.0;
}

static inline void
lb_d3q19_BGK_collide_pdfs_1(LBD3Q19BGK* self, double ctau, int o)
{
	double eq[19];
	const double ohm = 1./ctau;
	const double alp = 1. - 1./ctau;

	lb_d3q19_BGK_equilibrium(&self->parameters, self->rho + o,
				   self->u + 3*o, eq);

	{
		int n;

		for (n = 0; n < 19; ++n) {
			const double tmp = self->pdfs[n][o];
			self->pdfs[n][o] = alp*tmp + ohm*eq[n];
		}
	}
}

static void
lb_d3q19_BGK_collide_pdfs(LBD3Q19BGK* self, double ctau)
{
	static const double zero = 0.0;
	int x, y, z, z0, z1;

	lb_assert(self != NULL);

#ifdef LB_ENABLE_MPI
	z1 = self->d3q19.neighbors[4] == MPI_PROC_NULL ? (NZ - 1) : NZ;
	z0 = self->d3q19.neighbors[5] == MPI_PROC_NULL ? 2 : 1;
#else
	if (self->d3q19.partition_info.periods[2]) {
		z0 = 1;
		z1 = NZ;
	} else {
		z0 = 2;
		z1 = NZ - 1;
	}
#endif /* LB_ENABLE_MPI */

	/*
	 *  bulk sites first
	 */

	for (x = 1; x <= NX; ++x)
		for (y = 1; y <= NY; ++y)
			for (z = z0; z <= z1; ++z) {
				const int o = OFFSET(x, y, z);
				const int dox = (NY + 2)*(NZ + 2);
				const int doy = NZ + 2;

				const double drx = 0.5*(self->rho[o + dox]
						      - self->rho[o - dox]);
				const double dry = 0.5*(self->rho[o + doy]
						      - self->rho[o - doy]);
				const double drz = 0.5*(self->rho[o + 1]
						      - self->rho[o - 1]);
				const double ddr =  self->rho[o + dox]
						  + self->rho[o - dox]
						  + self->rho[o + doy]
						  + self->rho[o - doy]
						  + self->rho[o + 1]
						  + self->rho[o - 1]
						  - 6.0*self->rho[o];

				lb_d3q19_BGK_collide_pdfs_1(self, ctau, o);
			}

	/*
	 *  boundaries
	 */

	if (z0 == 2)
		for (x = 1; x <= NX; ++x)
			for (y = 1; y <= NY; ++y) {
				const int o = OFFSET(x, y, 1);

				const double drx =
				    0.5*(self->rho[OFFSET(x + 1, y, 1)]
				       - self->rho[OFFSET(x - 1, y, 1)]);

				const double dry =
				    0.5*(self->rho[OFFSET(x, y + 1, 1)]
				       - self->rho[OFFSET(x, y - 1, 1)]);

				const double ddr =
				    self->rho[OFFSET(x + 1, y, 1)]
				  + self->rho[OFFSET(x - 1, y, 1)]
				  + self->rho[OFFSET(x, y + 1, 1)]
				  + self->rho[OFFSET(x, y - 1, 1)]
				  + 0.5*self->rho[OFFSET(x, y, 2)]
				  - 4.5*self->rho[o];

				lb_d3q19_BGK_apply_bc(self, o, self->uxuy_b5,
							10, 11, 14, 15);

				lb_d3q19_BGK_collide_pdfs_1(self, ctau, o);
			}

	if (z1 == NZ - 1)
		for (x = 1; x <= NX; ++x)
			for (y = 1; y <= NY; ++y) {
				const int o = OFFSET(x, y, NZ);

				const double drx =
				    0.5*(self->rho[OFFSET(x + 1, y, NZ)]
				       - self->rho[OFFSET(x - 1, y, NZ)]);

				const double dry =
				    0.5*(self->rho[OFFSET(x, y + 1, NZ)]
				       - self->rho[OFFSET(x, y - 1, NZ)]);

				const double ddr =
				    self->rho[OFFSET(x + 1, y, NZ)]
				  + self->rho[OFFSET(x - 1, y, NZ)]
				  + self->rho[OFFSET(x, y + 1, NZ)]
				  + self->rho[OFFSET(x, y - 1, NZ)]
				  + 0.5*self->rho[OFFSET(x, y, NZ - 1)]
				  - 4.5*self->rho[o];

				lb_d3q19_BGK_apply_bc(self, o, self->uxuy_b4,
							12, 13, 14, 17);

				lb_d3q19_BGK_collide_pdfs_1(self, ctau, o);
			}
}

#ifdef LB_ENABLE_MPI
static void
lb_d3q19_BGK_exchange_rho(LBD3Q19BGK* self)
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
		if (self->d3q19.neighbors[n] != MPI_PROC_NULL &&
		    self->d3q19.neighbors[n] != THIS_PROCESS) {
			MPI_Datatype type = MPI_DATATYPE_NULL;

			switch (n/2) {
				case 0 : type = self->d3q19.plane_YZ_type; break;
				case 1 : type = self->d3q19.plane_XZ_type; break;
				case 2 : type = self->d3q19.plane_XY_type; break;
				default : lb_assert(FALSE);
			}

			MPI_Isend(self->rho + OFFSET_V(todo[n].s), 1, type,
				  self->d3q19.neighbors[n], n,
				  self->d3q19.comm, self->d3q19.requests + r++);

			MPI_Irecv(self->rho + OFFSET_V(todo[n].r), 1, type,
				  self->d3q19.neighbors[n], opposite[n],
				  self->d3q19.comm, self->d3q19.requests + r++);
		}

	lb_assert(r <= self->d3q19.n_requests);

	if (r > 0) {
		MPI_Status* statuses = lb_newa(MPI_Status, r);
		MPI_Waitall(r, self->d3q19.requests, statuses);
	}
}
#endif /* LB_ENABLE_MPI */

static inline void
lb_d3q19_BGK_copy_X_edge(LBD3Q19BGK* self, double* dst, const double* src)
{
	const int stride = (NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		*dst = *src;
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_BGK_copy_Y_edge(LBD3Q19BGK* self, double* dst, const double* src)
{
	const int stride = NZ + 2;
	int count = NY;

	for (; count > 0; --count) {
		*dst = *src;
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_BGK_copy_Z_edge(LBD3Q19BGK* self, double* dst, const double* src)
{
	memcpy(dst, src, NZ*sizeof(double));
}

static inline void
lb_d3q19_BGK_copy_XY_plane(LBD3Q19BGK* self, double* dst, const double* src)
{
	const int stride = (NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		lb_d3q19_BGK_copy_Y_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_BGK_copy_XZ_plane(LBD3Q19BGK* self, double* dst, const double* src)
{
	const int stride = (NY + 2)*(NZ + 2);
	int count = NX;

	for (; count > 0; --count) {
		lb_d3q19_BGK_copy_Z_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static inline void
lb_d3q19_BGK_copy_YZ_plane(LBD3Q19BGK* self, double* dst, const double* src)
{
	const int stride = NZ + 2;
	int count = NY;

	for (; count > 0; --count) {
		lb_d3q19_BGK_copy_Z_edge(self, dst, src);
		dst += stride;
		src += stride;
	}
}

static void
lb_d3q19_BGK_wrap_rho_X(LBD3Q19BGK* self)
{
	lb_d3q19_BGK_copy_YZ_plane(self, self->rho + OFFSET(0, 1, 1),
					   self->rho + OFFSET(NX, 1, 1));
	lb_d3q19_BGK_copy_YZ_plane(self, self->rho + OFFSET(NX + 1, 1, 1),
					   self->rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_BGK_wrap_rho_Y(LBD3Q19BGK* self)
{
	lb_d3q19_BGK_copy_XZ_plane(self, self->rho + OFFSET(1, 0, 1),
					   self->rho + OFFSET(1, NY, 1));
	lb_d3q19_BGK_copy_XZ_plane(self, self->rho + OFFSET(1, NY + 1, 1),
					   self->rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_BGK_wrap_rho_Z(LBD3Q19BGK* self)
{
	lb_d3q19_BGK_copy_XY_plane(self, self->rho + OFFSET(1, 1, 0),
					   self->rho + OFFSET(1, 1, NZ));
	lb_d3q19_BGK_copy_XY_plane(self, self->rho + OFFSET(1, 1, NZ + 1),
					   self->rho + OFFSET(1, 1, 1));
}

static void
lb_d3q19_BGK_wrap_rho(LBD3Q19BGK* self)
{
#ifdef LB_ENABLE_MPI
        if (self->d3q19.neighbors[0] == THIS_PROCESS)
#else
        if (self->d3q19.partition_info.periods[0])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_BGK_wrap_rho_X(self);

#ifdef LB_ENABLE_MPI
        if (self->d3q19.neighbors[2] == THIS_PROCESS)
#else
        if (self->d3q19.partition_info.periods[1])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_BGK_wrap_rho_Y(self);

#ifdef LB_ENABLE_MPI
        if (self->d3q19.neighbors[4] == THIS_PROCESS)
#else
        if (self->d3q19.partition_info.periods[2])
#endif /* LB_ENABLE_MPI */
		lb_d3q19_BGK_wrap_rho_Z(self);
}

static void
lb_d3q19_BGK_setup_rho_guards(LBD3Q19BGK* self)
{
#ifdef LB_ENABLE_MPI
	lb_d3q19_BGK_exchange_rho(self);
#endif /* LB_ENABLE_MPI */
	lb_d3q19_BGK_wrap_rho(self);
}

/*
 *  public interface
 */

LBD3Q19BGK*
LBD3Q19BGK_new(int nx, int ny, int nz, int pz)
{
	int nnodes, n;
	LBD3Q19BGK* self;

	lb_assert(nx > 2 && ny > 2 && nz > 2);

	self = lb_new(LBD3Q19BGK, 1);

#ifdef LB_ENABLE_MPI
	lb_d3q19_init(&self->d3q19, MPI_COMM_WORLD, 1,
			nx, ny, nz, TRUE, TRUE, pz);
#else
	lb_d3q19_init(&self->d3q19, 1, nx, ny, nz, TRUE, TRUE, pz);
#endif /* LB_ENABLE_MPI */

	self->parameters.T = 1.0;
	self->parameters.a = 0.0;
	self->parameters.b = 0.0;
	self->parameters.K = 0.01;
	self->parameters.tau = 0.8;

	if (NX < 3 || NY < 3 || NZ < 3)
		lb_fatal("D3Q19BGK: size of the grid too small"
#ifdef LB_ENABLE_MPI
			   " (or too much MPI processes)"
#endif /* LB_ENABLE_MPI */
			  );

	nnodes = (NX + 2)*(NY + 2)*(NZ + 2);

	for (n = 0; n < 19; ++n)
		self->pdfs[n] = lb_new(double, nnodes);

	self->rho = lb_new(double, nnodes);
	self->u = lb_new(double, 3*nnodes);

	self->uxuy_b4[0] = self->uxuy_b4[1] =
	self->uxuy_b5[0] = self->uxuy_b5[1] = 0.0;

#ifdef LB_ENABLE_MPI
	LB_CALL_MPI(MPI_Barrier(self->d3q19.comm))
#endif /* LB_ENABLE_MPI */

	return self;
}

void
LBD3Q19BGK_destroy(LBD3Q19BGK* self)
{
	int n;
	lb_assert(self != NULL);

	lb_free(self->u);
	lb_free(self->rho);

	for (n = 0; n < 19; ++n)
		lb_free(self->pdfs[n]);

	lb_d3q19_fini(&self->d3q19);

	lb_free(self);
}

void
_LBD3Q19BGK_partition_info(const LBD3Q19BGK* self,
			      LBPartitionInfo* pinfo)
{
	lb_assert(self != NULL);
	lb_assert(pinfo != NULL);

	*pinfo = self->d3q19.partition_info;
}

LBPartitionInfo*
LBD3Q19BGK_partition_info(const LBD3Q19BGK* self)
{
	LBPartitionInfo* pinfo;

	lb_assert(self != NULL);

	pinfo = lb_new(LBPartitionInfo, 1);

	_LBD3Q19BGK_partition_info(self, pinfo);

	return pinfo;
}

void
LBD3Q19BGK_set_parameters(LBD3Q19BGK*            self,
			      const LBGKParameters* parameters)
{
	lb_assert(self != NULL);
	lb_assert(parameters != NULL);
	lb_assert(parameters->tau > 0.5);

	self->parameters = *parameters;
}

void
_LBD3Q19BGK_get_parameters(const LBD3Q19BGK* self,
			      LBGKParameters*  parameters)
{
	lb_assert(self != NULL);
	lb_assert(parameters != NULL);

	*parameters = self->parameters;
}

LBGKParameters*
LBD3Q19BGK_get_parameters(const LBD3Q19BGK* self)
{
	LBGKParameters* params;

	lb_assert(self != NULL);

	params = lb_new(LBGKParameters, 1);

	_LBD3Q19BGK_get_parameters(self, params);

	return params;
}


void
LBD3Q19BGK_set_walls_speed(LBD3Q19BGK* self,
			       const double v4[2],
			       const double v5[2])
{
	lb_assert(self != NULL);
	lb_assert(v4 != NULL && v5 != NULL);

	if (self->d3q19.partition_info.periods[2])
		lb_warning("LBD3Q19BGK: setting speed of the walls will"
			     " not have any effect because of periodic boundary"
			     " conditions along Z axis\n");

	self->uxuy_b4[0] = v4[0];
	self->uxuy_b4[1] = v4[1];

	self->uxuy_b5[0] = v5[0];
	self->uxuy_b5[1] = v5[1];
}

void
LBD3Q19BGK_get_walls_speed(const LBD3Q19BGK* self,
			       double v4[2], double v5[2])
{
	lb_assert(self != NULL);
	lb_assert(v4 != NULL && v5 != NULL);

	v4[0] = self->uxuy_b4[0];
	v4[1] = self->uxuy_b4[1];

	v5[0] = self->uxuy_b5[0];
	v5[1] = self->uxuy_b5[1];
}

void 
LBD3Q19BGK_set_pdfs(LBD3Q19BGK* self, int a, int x, int y, int z, double f)
{
	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(z >= 0 && z < NZ);
	lb_assert(f >= 0.0);
	
	self->pdfs[a][OFFSET(x + 1, y + 1, z + 1)] = f;
}

void 
LBD3Q19BGK_get_pdfs(const LBD3Q19BGK* self, int a, int x, int y, int z, double* f)
{
	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(z >= 0 && z < NZ);
	lb_assert(f);
	
	*f = self->pdfs[a][OFFSET(x + 1, y + 1, z + 1)];
}

void
LBD3Q19BGK_set_averages(LBD3Q19BGK* self, int x, int y, int z,
			    double rho, double ux, double uy, double uz)
{
	int o;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(z >= 0 && z < NZ);
	lb_assert(rho > 0.0);

	o = OFFSET(x + 1, y + 1, z + 1);
	self->rho[o] = rho;

	o *= 3;
	self->u[o + 0] = ux;
	self->u[o + 1] = uy;
	self->u[o + 2] = uz;
}

void
LBD3Q19BGK_get_averages(const LBD3Q19BGK* self, int x, int y, int z,
			    double* rho, double* ux, double* uy, double* uz)
{
	int o;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < NX);
	lb_assert(y >= 0 && y < NY);
	lb_assert(z >= 0 && z < NZ);
	lb_assert(rho && ux && uy && uz);

	o = OFFSET(x + 1, y + 1, z + 1);
	*rho = self->rho[o];

	o *= 3;
	*ux  = self->u[o + 0];
	*uy  = self->u[o + 1];
	*uz  = self->u[o + 2];
}

void
LBD3Q19BGK_set_equilibrium(LBD3Q19BGK* self)
{
	int n, N;
	double* ptr;

	lb_assert(self != NULL);

	lb_d3q19_BGK_setup_rho_guards(self);

	N = (NY + 2)*(NX + 2)*(NZ + 2);
	ptr = self->pdfs[0];
	for (n = N; n > 0; --n)
		*ptr++ = 0.0;

	for (n = 1; n < 19; ++n)
		memcpy(self->pdfs[n], self->pdfs[0], N*sizeof(double));

	lb_d3q19_BGK_collide_pdfs(self, 1.0);	
}

void
LBD3Q19BGK_advance(LBD3Q19BGK* self)
{
	lb_assert(self != NULL);

	lb_d3q19_stream_pdfs(&self->d3q19, self->pdfs);
	lb_d3q19_BGK_average_pdfs(self);
	lb_d3q19_BGK_setup_rho_guards(self);
	lb_d3q19_BGK_collide_pdfs(self, self->parameters.tau);
}

void
LBD3Q19BGK_dump(const LBD3Q19BGK* self, const char* filename)
{
	hid_t file_id;
	hid_t file_space_id;
	hid_t mem_space_id;
	hid_t vel_space_id;

	hid_t rho_set_id;
	hid_t ux_set_id;
	hid_t uy_set_id;
	hid_t uz_set_id;

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

		H5Pset_fapl_mpio(acc_tpl, self->d3q19.comm, MPI_INFO_NULL);
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
		hsize_t gdims[3] = {self->d3q19.partition_info.global_size[0],
				    self->d3q19.partition_info.global_size[1],
				    self->d3q19.partition_info.global_size[2]};
#else
		hsize_t gdims[3] = {NX, NY, NZ};
#endif /* LB_ENABLE_MPI */

		hsize_t ldims[3] = {(NX + 2), (NY + 2), (NZ + 2)};
		hsize_t udims[3] = {(NX + 2), (NY + 2), 3*(NZ + 2)};

		file_space_id = H5Screate_simple(3, gdims, NULL);
		mem_space_id = H5Screate_simple(3, ldims, NULL);
		vel_space_id = H5Screate_simple(3, udims, NULL);
	}

	/*
	 *  select memory slab
	 */

	{
		hssize_t start[3] = {1, 1, 1};
		hsize_t  count[3] = {NX, NY, NZ};
		hsize_t  stride[3] = {1, 1, 1};

		hssize_t start2[3] = {1, 1, 3};
		hsize_t  count2[3] = {NX, NY, NZ};
		hsize_t  stride2[3] = {1, 1, 3};

		H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				    start, stride, count, NULL);
		H5Sselect_hyperslab(vel_space_id, H5S_SELECT_SET,
				    start2, stride2, count2, NULL);
	}

	/*
	 *  select slab in the file's data space
	 */

#ifdef LB_ENABLE_MPI
	{
		hssize_t start[3] = {self->d3q19.partition_info.global_origin[0],
				     self->d3q19.partition_info.global_origin[1],
				     self->d3q19.partition_info.global_origin[2]};
		hsize_t  count[3] = {NX, NY, NZ};

		H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET,
				    start, NULL, count, NULL);
	}
#endif /* LB_ENABLE_MPI */

	rho_set_id = H5Dcreate(file_id, "rho", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);
	ux_set_id = H5Dcreate(file_id, "ux", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);
	uy_set_id = H5Dcreate(file_id, "uy", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);
	uz_set_id = H5Dcreate(file_id, "uz", H5T_NATIVE_DOUBLE,
	                       file_space_id, H5P_DEFAULT);

#ifdef LB_ENABLE_MPI
	xfer_plist = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
#endif /* LB_ENABLE_MPI */

	H5Dwrite(rho_set_id, H5T_NATIVE_DOUBLE, mem_space_id,
		 file_space_id, xfer_plist, self->rho);
	H5Dwrite(ux_set_id, H5T_NATIVE_DOUBLE, vel_space_id,
		 file_space_id, xfer_plist, self->u);
	H5Dwrite(uy_set_id, H5T_NATIVE_DOUBLE, vel_space_id,
		 file_space_id, xfer_plist, self->u+1);
	H5Dwrite(uz_set_id, H5T_NATIVE_DOUBLE, vel_space_id,
		 file_space_id, xfer_plist, self->u+2);

#ifdef LB_ENABLE_MPI
	H5Pclose(xfer_plist);
#endif /* ICHF_ENABLE_MPI */

	H5Dclose(rho_set_id);
	H5Dclose(ux_set_id);
	H5Dclose(uy_set_id);
	H5Dclose(uz_set_id);

	H5Sclose(file_space_id);
	H5Sclose(mem_space_id);
	H5Sclose(vel_space_id);
	H5Fclose(file_id);
}

void
_LBD3Q19BGK_stats(const LBD3Q19BGK* self, LBD3Q19BGKStats* stats)
{
	int x, y, z;
	LBD3Q19BGKStats my;

	lb_assert(self != NULL);

	my.min_density = 1.0e15;
	my.max_density = my.max_velocity =
	my.kin_energy = my.momentum_x = my.momentum_y = my.momentum_z = 0.0;

	for (x = 1; x <= NX; ++x)
		for (y = 1; y <= NY; ++y)
			for (z = 1; z <= NZ; ++z) {
				int o = OFFSET(x, y, z);
				const double r  = self->rho[o];
				double ux, uy, uz, u2, u;

				o *= 3;
				ux = self->u[o++];
				uy = self->u[o++];
				uz = self->u[o];

				u2 = ux*ux + uy*uy;
				u = sqrt(u2);

				my.min_density = LB_MIN(my.min_density, r);
				my.max_density = LB_MAX(my.max_density, r);
				my.max_velocity = LB_MAX(my.max_velocity, u);
				my.kin_energy += r*u2;
				my.momentum_x += r*ux;
				my.momentum_y += r*uy;
				my.momentum_z += r*uz;
			}
#ifdef LB_ENABLE_MPI
        LB_CALL_MPI(MPI_Allreduce(&my.min_density, &stats->min_density,
		      1, MPI_DOUBLE, MPI_MIN, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.max_density, &stats->max_density,
		      1, MPI_DOUBLE, MPI_MAX, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.max_velocity,  &stats->max_velocity,
		      1, MPI_DOUBLE, MPI_MAX, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.kin_energy, &stats->kin_energy,
		      1, MPI_DOUBLE, MPI_SUM, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.momentum_x, &stats->momentum_x,
		      1, MPI_DOUBLE, MPI_SUM, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.momentum_y, &stats->momentum_y,
		      1, MPI_DOUBLE, MPI_SUM, self->d3q19.comm))
        LB_CALL_MPI(MPI_Allreduce(&my.momentum_z, &stats->momentum_z,
		      1, MPI_DOUBLE, MPI_SUM, self->d3q19.comm))

	{
		const int global_N =
			self->d3q19.partition_info.global_size[0]*
			self->d3q19.partition_info.global_size[1]*
			self->d3q19.partition_info.global_size[2];

		stats->kin_energy /= 2*global_N;
		stats->momentum_x /= global_N;
		stats->momentum_y /= global_N;
		stats->momentum_z /= global_N;
	}
#else /* !LB_ENABLE_MPI */
	my.kin_energy /= 2.0*NX*NY*NZ;
	my.momentum_x /= NX*NY*NZ;
	my.momentum_y /= NX*NY*NZ;
        *stats = my;
#endif /* LB_ENABLE_MPI */
}

LBD3Q19BGKStats*
LBD3Q19BGK_stats(const LBD3Q19BGK* self)
{
	LBD3Q19BGKStats* stat;

	lb_assert(self != NULL);

	stat = lb_new(LBD3Q19BGKStats, 1);

	_LBD3Q19BGK_stats(self, stat);

	return stat;
}

double
LBD3Q19BGK_mass(const LBD3Q19BGK* self)
{
	int x, y, z;
	double mymass = 0.0;

	lb_assert(self != NULL);

	for (x = 1; x <= NX; ++x)
		for (y = 1; y <= NY; ++y)
			for (z = 1; z <= NZ; ++z)
				mymass += self->rho[OFFSET(x, y, z)];

#ifdef LB_ENABLE_MPI
	{
		double totalmass;
		LB_CALL_MPI(MPI_Allreduce(&mymass, &totalmass, 1,
			      MPI_DOUBLE, MPI_SUM, self->d3q19.comm))
		return totalmass;
	}
#else
	return mymass;
#endif /* LB_ENABLE_MPI */
}
