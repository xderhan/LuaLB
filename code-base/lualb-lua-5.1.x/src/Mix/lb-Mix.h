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

#ifndef LB_Mix_H
#define LB_Mix_H

#include <lb-macros.h>

LB_BEGIN_DECLS

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

double lb_Mix_chemical_potential(const LBMixParameters*,
			    const double* rho,
			    const double* phi,
			    const double* ddp);

double lb_Mix_pressure_2d(const LBMixParameters*,
			    const double* rho,
				const double* phi,
			    const double* drx,
			    const double* dry,
				const double* dpx,
				const double* dpy,
				const double* ddr,
			    const double* ddp);

void lb_d2q9_Mix_equilibrium(const LBMixParameters*,
			       const double* rho,
				   const double* u,
			       const double* drx,
			       const double* dry,
				   const double* dpx,
			       const double* dpy,
				   const double* ddr,
			       const double* ddp,
			       double* equilibrium_pdfs[9]);

double lb_Mix_pressure_3d(const LBMixParameters*,
		     const double* rho,
			 const double* phi,
		     const double* drx,
		     const double* dry,
		     const double* drz,
			 const double* dpx,
		     const double* dpy,
		     const double* dpz,
			 const double* ddr,
		     const double* ddp);

void lb_d3q19_Mix_equilibrium(const LBMixParameters*,
			        const double* rho,
					const double* u,
					const double* drx,
					const double* dry,
					const double* drz,
					const double* dpx,
					const double* dpy,
					const double* dpz,
					const double* ddr,
					const double* ddp,
			        double* equilibrium_pdfs[19]);

LB_END_DECLS

#endif /* LB_Mix_H */
