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

%module lb
%include "typemaps.i"

%{

#define LB_LUA_WRAP_C

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef LB_ENABLE_MPI
#include "lb-mpi-utils.h"
#endif /* LB_ENABLE_MPI */

#include <stdlib.h>

#include "lb-d2q9-BGK.h"
#include "lb-d3q19-BGK.h"
#include "lb-d2q9-Mix.h"
#include "lb-d3q19-Mix.h"
#include "lb-d2q9-Landau.h"
#include "lb-d3q19-Landau.h"

#include "lb-wtime.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_RGB
#include "lb-rgb.h"
#endif /* LB_ENABLE_RGB */

/* You may do not need to include these header files
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <lb-messages.h>
#include <lb-BGK.h>
#include <lb-Mix.h>
#include <lb-Landau.h>
#include <lb-macros.h>
#include <lb-partition-info.h>
*/

/* extending for LBPartitionInfo array membrs */
int LBPartitionInfo_ndims(LBPartitionInfo *self){ return self->ndims;}
int LBPartitionInfo_processor_rank(LBPartitionInfo *self){ return self->processor_rank;}
int LBPartitionInfo_processor_size(LBPartitionInfo *self, int i){ return self->processors_size[i];}
int LBPartitionInfo_processor_coords(LBPartitionInfo *self, int i){ return self->processor_coords[i];}
int LBPartitionInfo_size(LBPartitionInfo *self, int i){ return self->size[i];}
int LBPartitionInfo_periods(LBPartitionInfo *self, int i){ return self->periods[i];}
int LBPartitionInfo_global_size(LBPartitionInfo *self, int i){ return self->global_size[i];}
int LBPartitionInfo_global_origin(LBPartitionInfo *self, int i){ return self->global_origin[i];}

%}

typedef struct {
	%immutable;
	/*
	int ndims;
	int processor_rank;
	//%apply int OUTPUT[ANY] {processors_size[3],processor_coords[3],size[3],periods[3],global_size[3],global_origin[3]};
	int processors_size[3];
	int processor_coords[3];

	%apply int INPUT[ANY] {int size[3]};
	int size[3];
	%clear int size[3];
	int periods[3];
	int global_size[3];
	int global_origin[3];
	//%clear processors_size[3],processor_coords[3],size[3],periods[3],global_size[3],global_origin[3];
	*/
} LBPartitionInfo;


%extend LBPartitionInfo {
	int ndims();
	int processor_rank();
	int processor_size(int i);
	int processor_coords(int i);
	int size(int i);
	int periods(int i);
	int global_size(int i);
	int global_origin(int i);
}

typedef struct {
	%mutable;
	double tau;
} LBGKParameters;

typedef struct {
	double T;
	double a;
	double b;
	double K;
	double Gr;
	double lamda;
	double rtau;
	double ptau;
} LBMixParameters;

%extend LBMixParameters {
	void set_T(double x);
	void set_a(double x);
	void set_b(double x);
	void set_K(double x);
	void set_Gr(double x);
	void set_lamda(double x);
	void set_rtau(double x);
	void set_ptau(double x);
	void set(double, double, double, double, double, double, double, double);
}

typedef struct {
	double a;
	double b;
	double K;
	double rtau;
	double ptau;
} LBLandauParameters;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
} LBD2Q9BGKStats;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
	double momentum_z;
} LBD3Q19BGKStats;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
} LBD2Q9MixStats;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
	double momentum_z;
} LBD3Q19MixStats;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
} LBD2Q9LandauStats;

typedef struct {
	%immutable;
	double min_density;
	double max_density;
	double max_velocity;
	double kin_energy;
	double momentum_x;
	double momentum_y;
	double momentum_z;
} LBD3Q19LandauStats;

%nodefaultctor;
typedef struct {
	%immutable;
} LBD2Q9BGK;

%nodefaultctor;
typedef struct {
	%immutable;
} LBD3Q19BGK;

%nodefaultctor;
typedef struct {
	%immutable;
} LBD2Q9Mix;

%nodefaultctor;
typedef struct {
	%immutable;
} LBD3Q19Mix;

%nodefaultctor;
typedef struct {
	%immutable;
} lb_d2q9_Landau;

%nodefaultctor;
typedef struct {
	%immutable;
} lb_d3q19_Landau;

#ifdef LB_ENABLE_RGB

%nodefaultctor;
typedef struct {
	%immutable;
} lb_colormap;

%nodefaultctor;
typedef struct {
	%immutable;
} lb_rgb;

#endif /* LB_ENABLE_RGB */

%rename(d2q9_BGK) LBD2Q9BGK_new;
LBD2Q9BGK* LBD2Q9BGK_new(int c, int nx, int ny, int px, int py);
%rename(d3q19_BGK) LBD3Q19BGK_new;
LBD3Q19BGK* LBD3Q19BGK_new(int nx, int ny, int nz, int pz);

%rename(d2q9_Mix) LBD2Q9Mix_new;
LBD2Q9Mix* LBD2Q9Mix_new(int nx, int ny, int py);
%rename(d3q19_Mix) LBD3Q19Mix_new;
LBD3Q19Mix* LBD3Q19Mix_new(int nx, int ny, int nz, int pz);

%rename(d2q9_LD) lb_d2q9_Landau_new;
lb_d2q9_Landau* lb_d2q9_Landau_new(int nx, int ny, int py);
%rename(d3q19_LD) lb_d3q19_Landau_new;
lb_d3q19_Landau* lb_d3q19_Landau_new(int nx, int ny, int nz, int pz);

#ifdef LB_ENABLE_RGB

%rename(colormap) lb_colormap_new;
lb_colormap* lb_colormap_new(void);

%rename(rgb) lb_rgb_new;
lb_rgb* lb_rgb_new(int width, int height);

#endif /* LB_ENABLE_RGB */

/* ===== BGK module ===== */

%extend LBD2Q9BGK {

	void destroy();
	void partition_info(LBPartitionInfo* INOUT);
	//LBPartitionInfo* partition_info();
	LBD2Q9BGKStats* stats();
	//%apply TYPE *INOUT {const LBGKParameters *par};
	void set_parameters(const LBGKParameters *INPUT);
	//%clear const LBGKParameters *par;
	
	//%apply SWIGTYPE *OUTPUT[ANY] {const LBGKParameters *OUTPUT};
	void get_parameters(LBGKParameters *OUTPUT);
	//%clear const LBGKParameters *OUTPUT;
		
	//LBGKParameters* get_parameters();
	
	void set_walls_speed(double top, double bottom);
	void get_walls_speed(double *OUTPUT, double *OUTPUT);
	
	void set_pdfs(int a, int x, int y, int c, double f);	
	void get_pdfs(int a, int x, int y, int c, double *OUTPUT);

	void set_averages(int x, int y, int c, double density, double vx, double vy);
	void get_averages(int x, int y, int c, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	void set_equilibrium();
	
	void collide_node(int, int, int, int, double);
	void stream();
	void average();

	void advance();

	void dump(const char* filename);
	
	void diff_opt(int, int, int, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	double mass();

}

%extend LBD3Q19BGK {

	void destroy();
	//void partition_info(LBPartitionInfo* OUTPUT);
	LBPartitionInfo* partition_info();
	LBD3Q19BGKStats* stats();
	void set_parameters(const LBGKParameters*);
	LBGKParameters* get_parameters();
	
	%apply double INPUT[ANY] {const double v4[2],const double v5[2]};
	void set_walls_speed(const double v4[2], const double v5[2]);
	%clear const double v4[2],const double v5[2];
	%apply double OUTPUT[ANY] {double OUTPUT[2],double OUTPUT[2]};
	void get_walls_speed(double OUTPUT[2], double OUTPUT[2]);
	%clear double OUTPUT[2],double OUTPUT[2];	
	
	//void set_walls_speed(const double, const double, const double, const double);
	//void get_walls_speed(double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	
	void set_pdfs(int a, int x, int y, int z, double f);	
	void get_pdfs(int a, int x, int y, int z, double *OUTPUT);
	void set_averages(int x, int y, int z, double density, double vx, double vy, double vz);
	void get_averages(int x, int y, int z, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	void set_equilibrium();
	void advance();
	void dump(const char* filename);
	double mass();
}

/* ===== Bibary Mixture module ===== */

%extend LBD2Q9Mix {

	void destroy();

	LBPartitionInfo* partition_info();

	LBD2Q9MixStats* stats();

	void set_parameters(const LBMixParameters*);

	LBMixParameters* get_parameters();

	void set_walls_speed(double top, double bottom);

	void get_walls_speed(double *OUTPUT, double *OUTPUT);

	void set_averages(int x, int y, double density, double diff_density, double vx, double vy);

	void get_averages(int x, int y, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	void set_equilibrium();

	void advance();

	void dump(const char* filename);

	double mass();

}

%extend LBD3Q19Mix {

	void destroy();
	LBPartitionInfo* partition_info();
	LBD3Q19MixStats* stats();
	void set_parameters(const LBMixParameters*);
	LBMixParameters* get_parameters();

	//void set_walls_speed(const double, const double, const double, const double);
	%apply double INPUT[ANY] {const double v4[2],const double v5[2]};
	void set_walls_speed(const double v4[2], const double v5[2]);
	%clear const double v4[2],const double v5[2];
	%apply (double OUTPUT[ANY],double OUTPUT[ANY]) {(const double OUTPUT[2],const double OUTPUT[2])};
	void get_walls_speed(double OUTPUT[2], double OUTPUT[2]);
	%clear (const double OUTPUT[2],const double OUTPUT[2]);	
	//void get_walls_speed(double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	void set_averages(int x, int y, int z, double density, double diff_density, double vx, double vy, double vz);
	void get_averages(int x, int y, int z, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	void set_equilibrium();
	void advance();
	void dump(const char* filename);
	double mass();
}

/* ===== Binary Mixture Landau module ===== */
%extend lb_d2q9_Landau {

	void destroy();
	LBPartitionInfo* partition_info();
	LBD2Q9LandauStats* stats();
	void set_parameters(const LBLandauParameters*);
	LBLandauParameters* get_parameters();
	void set_walls_speed(double top, double bottom);
	void get_walls_speed(double *OUTPUT, double *OUTPUT);
	void set_averages(int x, int y, double density, double diff_density, double vx, double vy);
	void get_averages(int x, int y, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	void set_equilibrium();
	void advance();
	void dump(const char* filename);
	double mass();
}

%extend lb_d3q19_Landau {

	void destroy();
	LBPartitionInfo* partition_info();
	LBD3Q19LandauStats* stats();
	void set_parameters(const LBLandauParameters*);
	LBLandauParameters* get_parameters();
	%apply double INPUT[ANY] {const double v4[2],const double v5[2]};
	void set_walls_speed(const double v4[2], const double v5[2]);
	%clear const double v4[2],const double v5[2];
	%apply double OUTPUT[ANY] {double OUTPUT[2],double OUTPUT[2]};
	void get_walls_speed(double OUTPUT[2], double OUTPUT[2]);
	%clear double OUTPUT[2],double OUTPUT[2];	
	void set_averages(int x, int y, int z, double density, double diff_density, double vx, double vy, double vz);
	void get_averages(int x, int y, int z, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	void set_equilibrium();
	void advance();
	void dump(const char* filename);
	double mass();
}

/* ===== Utilities ===== */
%rename(wtime) lb_wtime;
double lb_wtime(void);

%rename(wtime_string) lb_wtime_string;
const char* lb_wtime_string(double);

%rename(is_parallel) lb_is_parallel;
static int lb_is_parallel(void);

%rename(info_enable) lb_info_enable;
void lb_info_enable(void);

%rename(info_disable) lb_info_disable;
void lb_info_disable(void);


#ifdef LB_ENABLE_MPI

%rename(mpi_barrier) lb_mpi_barrier;
void lb_mpi_barrier(void);

#endif /* LB_ENABLE_MPI */

/* ===== Graphics ===== */
#ifdef LB_ENABLE_RGB

/* Color map */
%extend lb_colormap {
	
	void destroy();
	int num_colors();
	%apply double INPUT[ANY] {const double c[3]};
	void set_color(int, const double c[3]);	void append_color(const double c[3]);
	%clear const double c[3];
	%apply double OUTPUT[ANY] {double OUTPUT[3]};
	void get_color(int, double OUTPUT[3]);
	void map_value(double, double OUTPUT[3]);
	%clear double OUTPUT[3];
}

/* RGB */
%extend lb_rgb {

	void destroy();
	int width();
	int height();
	%apply double INPUT[ANY] {const double rgb[3]};
	void fill(const double rgb[3]);
	void set_pixel(int x, int y, const double rgb[3]);
	%clear const double rgb[3];
	%apply double INPUT[ANY] {const double rgba[4]};
	void set_pixel_rgba(int, int, const double rgba[4]);
	%clear const double rgba[4];
	%apply double OUTPUT[ANY] {const double OUTPUT[3]};
	void get_pixel(int x, int y, double OUTPUT[3]);
	%clear const double OUTPUT[3];
	void map_value(const lb_colormap*, int, int, double);
	void save(char*);
}

#endif /* LB_ENABLE_RGB */
