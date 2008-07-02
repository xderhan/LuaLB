#ifndef LB_D2Q9_BGK_H
#define LB_D2Q9_BGK_H

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
} LBD2Q9BGKStats;

typedef struct _LBD2Q9BGK LBD2Q9BGK;

LBD2Q9BGK* LBD2Q9BGK_new(int c, int nx, int ny, int px, int py);
void LBD2Q9BGK_destroy(LBD2Q9BGK*);

LBPartitionInfo* LBD2Q9BGK_partition_info(const LBD2Q9BGK*);

void LBD2Q9BGK_set_parameters(LBD2Q9BGK*, const LBGKParameters*);

LBGKParameters* LBD2Q9BGK_get_parameters(const LBD2Q9BGK*);

void LBD2Q9BGK_set_walls_speed(LBD2Q9BGK*, double top, double bottom);
void LBD2Q9BGK_get_walls_speed(const LBD2Q9BGK*, double*, double*);

void LBD2Q9BGK_set_pdfs(LBD2Q9BGK*, int a, int x, int y, int c, double f);
void LBD2Q9BGK_get_pdfs(const LBD2Q9BGK*, int a, int x, int y, int c, double* f);

void LBD2Q9BGK_set_averages(LBD2Q9BGK*, int x, int y, int c, double density,
				double velocity_x, double velocity_y);
void LBD2Q9BGK_get_averages(const LBD2Q9BGK*, int x, int y, int c,
				double* density,
				double* velocity_x, double* velocity_y);

void LBD2Q9BGK_set_equilibrium(LBD2Q9BGK*);
void LBD2Q9BGK_collide_node(LBD2Q9BGK*, int, int, int, int, double);
void LBD2Q9BGK_stream(LBD2Q9BGK*);
void LBD2Q9BGK_average(LBD2Q9BGK*);
void LBD2Q9BGK_advance(LBD2Q9BGK*);
void LBD2Q9BGK_dump(const LBD2Q9BGK*, const char* filename);

LBD2Q9BGKStats* LBD2Q9BGK_stats(const LBD2Q9BGK*);

// remove this
double LBD2Q9BGK_mass(const LBD2Q9BGK*);
void LBD2Q9BGK_diff_opt(LBD2Q9BGK*, int, int, int, double *, double *, double *);

LB_END_DECLS

#endif /* LB_D2Q9_BGK_H */
