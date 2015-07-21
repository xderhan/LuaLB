#ifndef LB_D3Q19_H
#define LB_D3Q19_H

/*
 *  config.h should be already included in order to have proper LB_ENABLE_MPI
 */

#ifdef LB_ENABLE_MPI
#include <mpi.h>
#endif /* LB_ENABLE_MPI */

#include <lb-macros.h>
#include <lb-partition-info.h>

LB_BEGIN_DECLS

/*
 *  velocities:
 *
 *  (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
 *
 *  (1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0),
 *  (1, 0, 1), (-1, 0, 1), (1, 0, -1), (-1, 0, -1),
 *  (0, 1, 1), (0, -1, 1), (0, 1, -1), (0, -1, -1)
 *
 *  memory layout corresponds to that expected by HDF5 - last coordinate (Z)
 *  is fastest, i.e. PDFs (each of size particle_rank*sizeof(double)) are
 *  the in following order
 *
 *  (0, 0, 0), (0, 0, 1), (0, 0, 2) ... (0, 0, NZ + 1),
 *  (0, 1, 0), (0, 1, 1) ... and so on
 *
 */

typedef struct {
	LBPartitionInfo partition_info;
	int               particle_rank;
#ifdef LB_ENABLE_MPI
	MPI_Comm     comm;
	MPI_Datatype edge_X_type;
	MPI_Datatype edge_Y_type;
	MPI_Datatype edge_Z_type;
	MPI_Datatype plane_XY_type;
	MPI_Datatype plane_XZ_type;
	MPI_Datatype plane_YZ_type;
	int          neighbors[18];
	int          n_requests;
	MPI_Request* requests;
#endif /* LB_ENABLE_MPI */
} LBD3Q19;

void lb_d3q19_init(LBD3Q19*,
#ifdef LB_ENABLE_MPI
		     MPI_Comm comm,
#endif /* LB_ENABLE_MPI */
		     int prank,
		     int nx, int ny, int nz,
		     int px, int py, int pz);

void lb_d3q19_fini(LBD3Q19*);

void lb_d3q19_stream_pdfs(const LBD3Q19*, double* pdfs[18]);

LB_END_DECLS

#endif /* LB_D3Q19_H */
