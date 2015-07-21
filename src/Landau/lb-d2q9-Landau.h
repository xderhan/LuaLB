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

#ifndef LB_D2Q9_Landau_H
#define LB_D2Q9_Landau_H

#include "lb-Landau.h"

#include <lb-macros.h>
#include <lb-partition-info.h>

LB_BEGIN_DECLS

typedef struct {
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
} LBD2Q9LandauStats;

typedef struct _LBD2Q9Landau lb_d2q9_Landau;

lb_d2q9_Landau* lb_d2q9_Landau_new(int nx, int ny, int py);
void lb_d2q9_Landau_destroy(lb_d2q9_Landau*);

LBPartitionInfo* lb_d2q9_Landau_partition_info(const lb_d2q9_Landau*);
//void LBD2Q9VdW_partition_info(const LBD2Q9VdW*, LBPartitionInfo*);

void lb_d2q9_Landau_set_parameters(lb_d2q9_Landau*, const LBLandauParameters*);
//void LBD2Q9VdW_get_parameters(const LBD2Q9VdW*, LBVdWParameters*);
LBLandauParameters* lb_d2q9_Landau_get_parameters(const lb_d2q9_Landau*);

void lb_d2q9_Landau_set_walls_speed(lb_d2q9_Landau*, double top, double bottom);
void lb_d2q9_Landau_get_walls_speed(const lb_d2q9_Landau*, double*, double*);

void lb_d2q9_Landau_set_averages(lb_d2q9_Landau*, int x, int y, double density,
				double diff_density, double velocity_x, double velocity_y);
void lb_d2q9_Landau_get_averages(const lb_d2q9_Landau*, int x, int y,
				double* density, double* diff_density,
				double* velocity_x, double* velocity_y);

void lb_d2q9_Landau_set_equilibrium(lb_d2q9_Landau*);
void lb_d2q9_Landau_advance(lb_d2q9_Landau*);
void lb_d2q9_Landau_dump(const lb_d2q9_Landau*, const char* filename);

//void LBD2Q9VdW_stats(const LBD2Q9VdW*, LBD2Q9VdWStats*);
LBD2Q9LandauStats* lb_d2q9_Landau_stats(const lb_d2q9_Landau*);

// remove this
double lb_d2q9_Landau_mass(const lb_d2q9_Landau*);

LB_END_DECLS

#endif /* LB_D2Q9_Landau_H */
