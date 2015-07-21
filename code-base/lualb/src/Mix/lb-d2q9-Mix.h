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

#ifndef LB_D2Q9_Mix_H
#define LB_D2Q9_Mix_H

#include "lb-Mix.h"

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
} LBD2Q9MixStats;

typedef struct _LBD2Q9Mix LBD2Q9Mix;

LBD2Q9Mix* LBD2Q9Mix_new(int nx, int ny, int py);
void LBD2Q9Mix_destroy(LBD2Q9Mix*);

LBPartitionInfo* LBD2Q9Mix_partition_info(const LBD2Q9Mix*);
//void LBD2Q9VdW_partition_info(const LBD2Q9VdW*, LBPartitionInfo*);

void LBD2Q9Mix_set_parameters(LBD2Q9Mix*, const LBMixParameters*);
//void LBD2Q9VdW_get_parameters(const LBD2Q9VdW*, LBVdWParameters*);
LBMixParameters* LBD2Q9Mix_get_parameters(const LBD2Q9Mix*);

void LBD2Q9Mix_set_walls_speed(LBD2Q9Mix*, double top, double bottom);
void LBD2Q9Mix_get_walls_speed(const LBD2Q9Mix*, double*, double*);

void LBD2Q9Mix_set_averages(LBD2Q9Mix*, int x, int y, double density,
				double diff_density, double velocity_x, double velocity_y);
void LBD2Q9Mix_get_averages(const LBD2Q9Mix*, int x, int y,
				double* density, double* diff_density,
				double* velocity_x, double* velocity_y);

void LBD2Q9Mix_set_equilibrium(LBD2Q9Mix*);
void LBD2Q9Mix_advance(LBD2Q9Mix*);
void LBD2Q9Mix_dump(const LBD2Q9Mix*, const char* filename);

//void LBD2Q9VdW_stats(const LBD2Q9VdW*, LBD2Q9VdWStats*);
LBD2Q9MixStats* LBD2Q9Mix_stats(const LBD2Q9Mix*);

// remove this
double LBD2Q9Mix_mass(const LBD2Q9Mix*);

LB_END_DECLS

#endif /* LB_D2Q9_Mix_H */
