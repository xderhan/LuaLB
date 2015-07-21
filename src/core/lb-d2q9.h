#ifndef LB_D2Q9_H
#define LB_D2Q9_H

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
 *  velocities: v0 = (1, 0), v1 = (0, 1), v2 = (-1, 0), v3 = (0, -1),
 *              v4 = (1, 1), v5 = (-1, 1), v6 = (-1, -1), v7 = (-1, 1),
 *              v8 = (0, 0)
 *
 *  (the first 8 velocities coincide with the shifts to the MPI neighbors)
 *
 *  memory layout of PDFs is as follows
 *
 *  5   ^ Y                   MPI neighbor 1               MPI neighbor 4
 *    \ |                          ^                            /
 *      +----------------------- ..... -----------------------+
 *      |                   guard layer                       |
 *      |     +----------------- ..... -------------------+   |
 *      |     |                                   (NX, NY)|   |
 *      :     :                                           :   : --> 0
 * 2 <  :     :                                           :   :
 *      |     |(1,1)                                      |   |
 *      |     +----------------- ..... -------------------+   |
 *      |(0,0)|(0,1)|                                         |
 *      +------------------------- ..... ---------------------+ ---> X
 *    /                             |                          \
 *  neighbor 6                  MPI neighbor 3              MPI neighbor 7
 *
 *
 *  memory layout corresponds to that expected by HDF5 - last coordinate (Y)
 *  is fastest, i.e. PDFs (each of size particle_rank*sizeof(double)) are
 *  the in following order
 *
 *  (0, 0), (0, 1), (0, 2) ... (0, NY + 1), (1, 0), (1, 1) ... and so on
 *
 */

typedef struct {
	LBPartitionInfo partition_info;
	int             particle_rank;
#ifdef LB_ENABLE_MPI
	MPI_Comm     comm;
	MPI_Datatype col_type;
	MPI_Datatype row_type;
	int          neighbors[8];
	int          n_requests;
	MPI_Request* requests;
#endif /* LB_ENABLE_MPI */
} LBD2Q9;

void lb_d2q9_init(LBD2Q9*,
#ifdef LB_ENABLE_MPI
		    MPI_Comm comm,
#endif /* LB_ENABLE_MPI */
		    int prank,
		    int nx, int ny,
		    int px, int py);

void lb_d2q9_fini(LBD2Q9*);

void lb_d2q9_stream_pdfs(const LBD2Q9*, double* pdfs[8]);

#ifdef LB_ENABLE_MPI
/*
 *  returns either MPI_PROC_NULL or rank of the neighbor
 */
int lb_d2q9_need_exchange(const LBD2Q9*, int);
#endif /* LB_ENABLE_MPI */

int lb_d2q9_is_boundary(const LBD2Q9*, int);

LB_END_DECLS

#endif /* LB_D2Q9_H */
