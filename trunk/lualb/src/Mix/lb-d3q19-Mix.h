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

#ifndef LB_D3Q19_Mix_H
#define LB_D3Q19_Mix_H

#include <lb-Mix.h>
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
	double momentum_z;
} LBD3Q19MixStats;

typedef struct _LBD3Q19Mix LBD3Q19Mix;

LBD3Q19Mix* LBD3Q19Mix_new(int nx, int ny, int nz, int pz);
void LBD3Q19Mix_destroy(LBD3Q19Mix*);

//void LBD3Q19VdW_partition_info(const LBD3Q19VdW*, LBPartitionInfo*);
LBPartitionInfo* LBD3Q19Mix_partition_info(const LBD3Q19Mix*);

void LBD3Q19Mix_set_parameters(LBD3Q19Mix*, const LBMixParameters*);
//void LBD3Q19VdW_get_parameters(const LBD3Q19VdW*, LBVdWParameters*);
LBMixParameters* LBD3Q19Mix_get_parameters(const LBD3Q19Mix*);

void LBD3Q19Mix_set_walls_speed(LBD3Q19Mix*, const double v4[2],
						   const double v5[2]);
void LBD3Q19Mix_get_walls_speed(const LBD3Q19Mix*, double*, double*);

void LBD3Q19Mix_set_averages(LBD3Q19Mix*, int x, int y, int z,
    double density, double diff_density, double velocity_x, double velocity_y, double velocity_z);

void LBD3Q19Mix_get_averages(const LBD3Q19Mix*, int x, int y, int z,
    double* density, double* diff_density, double* velocity_x, double* velocity_y, double* velocity_z);

void LBD3Q19Mix_set_equilibrium(LBD3Q19Mix*);
void LBD3Q19Mix_advance(LBD3Q19Mix*);
void LBD3Q19Mix_dump(const LBD3Q19Mix*, const char* filename);

//void LBD3Q19VdW_stats(const LBD3Q19VdW*, LBD3Q19VdWStats*);
LBD3Q19MixStats* LBD3Q19Mix_stats(const LBD3Q19Mix*);

// remove this
double LBD3Q19Mix_mass(const LBD3Q19Mix*);

LB_END_DECLS

#endif /* LB_D3Q19_Mix_H */
