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

#ifndef LB_Landau_H
#define LB_Landau_H

#include <lb-macros.h>

LB_BEGIN_DECLS

typedef struct {
	double a;
	double b;
	double K;
	double rtau;
	double ptau;
} LBLandauParameters;

double lb_Landau_chemical_potential(const LBLandauParameters*,
			    const double* rho,
			    const double* phi,
			    const double* ddp);

double lb_Landau_pressure_2d(const LBLandauParameters*,
			    const double* rho,
				const double* phi,
			    const double* drx,
			    const double* dry,
				const double* dpx,
				const double* dpy,
				const double* ddr,
			    const double* ddp);

void lb_d2q9_Landau_equilibrium(const LBLandauParameters*,
			       const double* rho,
				   const double* u,
			       const double* drx,
			       const double* dry,
				   const double* dpx,
			       const double* dpy,
				   const double* ddr,
			       const double* ddp,
			       double* equilibrium_pdfs[9]);

double lb_Landau_pressure_3d(const LBLandauParameters*,
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

void lb_d3q19_Landau_equilibrium(const LBLandauParameters*,
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

#endif /* LB_Landau_H */
