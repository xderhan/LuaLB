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

#ifndef LB_D3Q19_Landau_H
#define LB_D3Q19_Landau_H

#include <lb-Landau.h>
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
} LBD3Q19LandauStats;

typedef struct _LBD3Q19Landau lb_d3q19_Landau;

lb_d3q19_Landau* lb_d3q19_Landau_new(int nx, int ny, int nz, int pz);
void lb_d3q19_Landau_destroy(lb_d3q19_Landau*);

LBPartitionInfo* lb_d3q19_Landau_partition_info(const lb_d3q19_Landau*);

void lb_d3q19_Landau_set_parameters(lb_d3q19_Landau*, const LBLandauParameters*);
LBLandauParameters* lb_d3q19_Landau_get_parameters(const lb_d3q19_Landau*);

void lb_d3q19_Landau_set_walls_speed(lb_d3q19_Landau*, const double v4[2],
						   const double v5[2]);
void lb_d3q19_Landau_get_walls_speed(const lb_d3q19_Landau*, double*, double*);

void lb_d3q19_Landau_set_averages(lb_d3q19_Landau*, int x, int y, int z,
    double density, double diff_density, double velocity_x, double velocity_y, double velocity_z);

void lb_d3q19_Landau_get_averages(const lb_d3q19_Landau*, int x, int y, int z,
    double* density, double* diff_density, double* velocity_x, double* velocity_y, double* velocity_z);

void lb_d3q19_Landau_set_equilibrium(lb_d3q19_Landau*);
void lb_d3q19_Landau_advance(lb_d3q19_Landau*);
void lb_d3q19_Landau_dump(const lb_d3q19_Landau*, const char* filename);

LBD3Q19LandauStats* lb_d3q19_Landau_stats(const lb_d3q19_Landau*);

// remove this
double lb_d3q19_Landau_mass(const lb_d3q19_Landau*);

LB_END_DECLS

#endif /* LB_D3Q19_Landau_H */
