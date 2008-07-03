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

#include "lb-wtime.h"
#include "lb-messages.h"

#ifdef LB_ENABLE_RGB
#include "lb-rgb.h"
#endif /* LB_ENABLE_RGB */

/* already added
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <lb-messages.h>
#include <lb-BGK.h>
#include <lb-Mix.h>
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

/* extending for LBGKParameters array membrs */
void LBGKParameters_set_tau(LBGKParameters *self, double x) {self->tau = x;}

void LBGKParameters_set(LBGKParameters *self, double tau) {
	self->tau = tau;
}

/* extending for LBMixParameters array membrs */
void LBMixParameters_set_T(LBMixParameters *self, double x) {self->T = x;}
void LBMixParameters_set_a(LBMixParameters *self, double x) {self->a = x;}
void LBMixParameters_set_b(LBMixParameters *self, double x) {self->b = x;}
void LBMixParameters_set_K(LBMixParameters *self, double x) {self->K = x;}
void LBMixParameters_set_Gr(LBMixParameters *self, double x) {self->Gr = x;}
void LBMixParameters_set_lamda(LBMixParameters *self, double x) {self->lamda = x;}
void LBMixParameters_set_rtau(LBMixParameters *self, double x) {self->rtau = x;}
void LBMixParameters_set_ptau(LBMixParameters *self, double x) {self->ptau = x;}

void LBMixParameters_set(LBMixParameters *self, 
				double T, 
				double a, 
				double b, 
				double K, 
				double Gr,
				double lamda,
				double rtau,
				double ptau) 
{
	self->T = T;
	self->a = a;
	self->b = b;
	self->K = K;
	self->Gr = Gr;
	self->lamda = lamda;
	self->rtau = rtau;
	self->ptau = ptau;
}

void LBD3Q19BGK_set_walls_speed(LBD3Q19BGK *self, const double v40, const double v41, 
						const double v50, const double v51) 
{
	double v4[2], v5[2];
	v4[0]=v41; v4[1]=v41; v5[0]=v50; v5[1]=v51;
	_LBD3Q19BGK_set_walls_speed(self, v4, v5);
}

void LBD3Q19BGK_get_walls_speed(LBD3Q19BGK *self, double *v40, double *v41, 
						double *v50, double *v51) 
{
	double v4[2], v5[2];
	_LBD3Q19BGK_get_walls_speed(self, v4, v5);
	*v40=v4[0]; *v41=v4[1]; *v50=v5[0]; *v51=v5[1];
}

void LBD3Q19Mix_set_walls_speed(LBD3Q19Mix *self, const double v40, const double v41, 
						const double v50, const double v51) 
{
	double v4[2], v5[2];
	v4[0]=v41; v4[1]=v41; v5[0]=v50; v5[1]=v51;
	_LBD3Q19Mix_set_walls_speed(self, v4, v5);
}

void LBD3Q19Mix_get_walls_speed(LBD3Q19Mix *self, double *v40, double *v41, 
						double *v50, double *v51) 
{
	double v4[2], v5[2];
	_LBD3Q19Mix_get_walls_speed(self, v4, v5);
	*v40=v4[0]; *v41=v4[1]; *v50=v5[0]; *v51=v5[1];
}


/* ========== Utilities ========== */

/* extending for LBColormap */
//void set_color(int, const double c[3]);
void LBColormap_set_color(const LBColormap* self, int c, 
				const double r, 
				const double g, 
				const double b) 
{
	double cl[3];
	cl[0]=r; cl[1]=g; cl[2]=b;
	_LBColormap_set_color(self, c, cl);
}

//void get_color(int, double c[3]);
void LBColormap_get_color(const LBColormap* self, int c, 
					double *r, 
					double *g, 
					double *b) 
{
	double cl[3];
	_LBColormap_get_color(self, c, cl);
	*r=cl[0]; *b=cl[1]; *b=cl[2];
}

//void append_color(const double c[3]);
/*
void LBColormap_append_color(LBColormap* self, 
				const double r, 
				const double g, 
				const double b) 
{
	double cl[3];
	cl[0]=r; cl[1]=g; cl[2]=b;
	_LBColormap_append_color(self, cl);
}
*/

//void map_value(double, double c[3]);
void LBColormap_map_value(const LBColormap* self, double v, 
				double *r, 
				double *g, 
				double *b) 
{
	double cl[3];
	_LBColormap_map_value(self, v, cl);
	*r=cl[0]; *b=cl[1]; *b=cl[2];
}

/* extending for LBRGB */
//void fill(const double rgb[3]);
void LBRGB_fill(LBRGB* self, 
			const double r, 
			const double g, 
			const double b) 
{
	double cl[3];
	cl[0]=r; cl[1]=g; cl[2]=b;
	_LBRGB_fill(self, cl);
}

//void get_pixel(int x, int y, double rgb[3]);
void LBRGB_get_pixel(const LBRGB* self, int x, int y, 
			double *r, 
			double *g, 
			double *b) 
{
	double cl[3];
	_LBRGB_get_pixel(self, x, y, cl);
	*r=cl[0]; *b=cl[1]; *b=cl[2];
}

//void set_pixel(int x, int y, const double rgb[3]);
void LBRGB_set_pixel(LBRGB* self, int x, int y, 
			const double r, 
			const double g, 
			const double b) 
{
	double cl[3];
	cl[0]=r; cl[1]=g; cl[2]=b;
	_LBRGB_set_pixel(self, x, y, cl);
}

//void set_pixel_rgba(int, int, const double rgba[4]);
void LBRGB_set_pixel_rgba(LBRGB* self, int x, int y, 
				const double r, 
				const double g, 
				const double b, 
				const double a) 
{
	double cl[4];
	cl[0]=r; cl[1]=g; cl[2]=b; cl[3]=a;
	_LBRGB_set_pixel_rgba(self, x, y, cl);
}

void LBRGB_save(LBRGB* self, char* name) {
	FILE *file;
	file = fopen(name, "w");
	_LBRGB_save(self, file);
}

%}

typedef struct {
	%immutable;
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
	double tau;
} LBGKParameters;

%extend LBGKParameters {
	void set_tau(double);
	void set(double);
}

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


#ifdef LB_ENABLE_RGB

%nodefaultctor;
typedef struct {
	%immutable;
} LBColormap;

%nodefaultctor;
typedef struct {
	%immutable;
} LBRGB;

#endif /* LB_ENABLE_RGB */

%rename(d2q9_BGK) LBD2Q9BGK_new;
LBD2Q9BGK* LBD2Q9BGK_new(int c, int nx, int ny, int px, int py);

%rename(d3q19_BGK) LBD3Q19BGK_new;
LBD3Q19BGK* LBD3Q19BGK_new(int nx, int ny, int nz, int pz);

%rename(d2q9_Mix) LBD2Q9Mix_new;
LBD2Q9Mix* LBD2Q9Mix_new(int nx, int ny, int py);

%rename(d3q19_Mix) LBD3Q19Mix_new;
LBD3Q19Mix* LBD3Q19Mix_new(int nx, int ny, int nz, int pz);


#ifdef LB_ENABLE_RGB

%rename(colormap) LBColormap_new;
LBColormap* LBColormap_new(void);

%rename(rgb) LBRGB_new;
LBRGB* LBRGB_new(int width, int height);

#endif /* LB_ENABLE_RGB */

/* ===== BGK module ===== */

%extend LBD2Q9BGK {

	void destroy();

	LBPartitionInfo* partition_info();

	LBD2Q9BGKStats* stats();

	void set_parameters(const LBGKParameters*);

	LBGKParameters* get_parameters();

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

	LBPartitionInfo* partition_info();

	LBD3Q19BGKStats* stats();

	void set_parameters(const LBGKParameters*);

	LBGKParameters* get_parameters();

	void set_walls_speed(const double, const double, const double, const double);

	void get_walls_speed(double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);
	
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

	void set_walls_speed(const double, const double, const double, const double);

	void get_walls_speed(double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);

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

%extend LBColormap {

	void destroy();

	int num_colors();

	//void set_color(int, const double c[3]);
	void set_color(int, const double, const double, const double);

	//void get_color(int, double c[3]);
	void get_color(int, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	%apply double INPUT[ANY] {const double c[3]};
	void append_color(const double c[3]);
	//void append_color(const double, const double, const double);
	%clear const double c[3];

	//void map_value(double, double c[3]);
	void map_value(double, double *OUTPUT, double *OUTPUT, double *OUTPUT);

}

/* RGB */

%extend LBRGB {

	void destroy();

	int width();

	int height();

	//void fill(const double rgb[3]);
	void fill(const double, const double, const double);

	//void get_pixel(int x, int y, double rgb[3]);
	void get_pixel(int x, int y, double *OUTPUT, double *OUTPUT, double *OUTPUT);

	//void set_pixel(int x, int y, const double rgb[3]);
	void set_pixel(int x, int y, const double, const double, const double);

	//void set_pixel_rgba(int, int, const double rgba[4]);
	void set_pixel_rgba(int, int, const double, const double, const double, const double);

	//void map_value(const LBColormap*, int, int, double);

	//void save(FILE*);
	void save(char*);

}

#endif /* LB_ENABLE_RGB */
