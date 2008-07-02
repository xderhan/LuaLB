#ifndef LB_D3Q19_BGK_H
#define LB_D3Q19_BGK_H

#include "lb-BGK.h"
#include "lb-macros.h"
#include "lb-partition-info.h"

LB_BEGIN_DECLS

typedef struct {
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
	double momentum_z;
} LBD3Q19BGKStats;

typedef struct _LBD3Q19BGK LBD3Q19BGK;

LBD3Q19BGK* LBD3Q19BGK_new(int nx, int ny, int nz, int pz);
void LBD3Q19BGK_destroy(LBD3Q19BGK*);

LBPartitionInfo* LBD3Q19BGK_partition_info(const LBD3Q19BGK*);

void LBD3Q19BGK_set_parameters(LBD3Q19BGK*, const LBGKParameters*);

LBGKParameters* LBD3Q19BGK_get_parameters(const LBD3Q19BGK*);

void _LBD3Q19BGK_set_walls_speed(LBD3Q19BGK*, const double v4[2],
						   const double v5[2]);
void _LBD3Q19BGK_get_walls_speed(const LBD3Q19BGK*, double*, double*);

void LBD3Q19BGK_set_pdfs(LBD3Q19BGK*, int a, int x, int y, int z, double f);
void LBD3Q19BGK_get_pdfs(const LBD3Q19BGK*, int a, int x, int y, int z, double* f);

void LBD3Q19BGK_set_averages(LBD3Q19BGK*, int x, int y, int z,
    double density, double velocity_x, double velocity_y, double velocity_z);

void LBD3Q19BGK_get_averages(const LBD3Q19BGK*, int x, int y, int z,
    double* density, double* velocity_x, double* velocity_y, double* velocity_z);

void LBD3Q19BGK_set_equilibrium(LBD3Q19BGK*);
void LBD3Q19BGK_advance(LBD3Q19BGK*);
void LBD3Q19BGK_dump(const LBD3Q19BGK*, const char* filename);

LBD3Q19BGKStats* LBD3Q19BGK_stats(const LBD3Q19BGK*);

// remove this
double LBD3Q19BGK_mass(const LBD3Q19BGK*);

LB_END_DECLS

#endif /* LB_D3Q19_BGK_H */
