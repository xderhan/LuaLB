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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <math.h>

#include "lb-Landau.h"
#include "lb-messages.h"

double lb_Landau_chemical_potential(const LBLandauParameters* params,
			    const double* LB_RESTRICT rho,
			    const double* LB_RESTRICT phi,
			    const double* LB_RESTRICT ddp)
{
	lb_assert(params && rho && phi && ddp);

	return params->a*(*phi) + params->b*(*phi)*(*phi)*(*phi) - params->K*(*ddp);
}

double
lb_Landau_pressure_2d(const LBLandauParameters* params,
		     const double* LB_RESTRICT rho,
			 const double* LB_RESTRICT phi,
		     const double* LB_RESTRICT drx,
		     const double* LB_RESTRICT dry,
			 const double* LB_RESTRICT dpx,
		     const double* LB_RESTRICT dpy,
		     const double* LB_RESTRICT ddr,
			 const double* LB_RESTRICT ddp)
{
	lb_assert(params && rho && phi && drx && dry && dpx && dpy && ddr && ddp);

	return (1./3.)*(*rho) + 0.5*params->a*(*phi)*(*phi) + 0.75*params->b*(*phi)*(*phi)*(*phi)*(*phi) 
			- params->K*( (*phi)*(*ddp) + 0.5*( (*dpx)*(*dpx) + (*dpy)*(*dpy) ) );
}

void
lb_d2q9_Landau_equilibrium(const LBLandauParameters* params,
			  const double* LB_RESTRICT averages,
			  const double* LB_RESTRICT u,
			  const double* LB_RESTRICT drx,
			  const double* LB_RESTRICT dry,
			  const double* LB_RESTRICT dpx,
			  const double* LB_RESTRICT dpy,
			  const double* LB_RESTRICT ddr,
			  const double* LB_RESTRICT ddp,
			  double** LB_RESTRICT eq)
{
	const double r = *averages++;
	const double p = *averages;
	const double ux = *u++;
	const double uy = *u;
	
	const double mu = lb_Landau_chemical_potential(params,&r,&p,ddp);
	const double Pr = lb_Landau_pressure_2d(params,&r,&p,drx,dry,dpx,dpy,ddr,ddp);
	
	const double muu = r*mu;
	const double Pxx = Pr + params->K*(*dpx)*(*dpx);
	const double Pyy = Pr + params->K*(*dpy)*(*dpy);
	const double Pxy =      params->K*(*dpx)*(*dpy);
	
	const double rux = r*ux;
	const double ruy = r*uy;
	const double ruxx = rux*ux + Pxx;
	const double ruyy = ruy*uy + Pyy;	
	const double ruxy = 2.*rux*uy + 2.*Pxy;	
	const double rusq = ruxx + ruyy;

	/* phi */
	const double pux = p*ux;
	const double puy = p*uy;
	const double puxx = pux*ux + muu;
	const double puyy = puy*uy + muu;	
	const double puxy = 2.*pux*uy;	
	const double pusq = puxx + puyy;

	/* rho eq dist. */
	/* v0 - v3 */
	{
		const double w = 1./9.;
		eq[0][0] = w*( 3.*rux+4.5*ruxx-1.5*rusq);
		eq[1][0] = w*( 3.*ruy+4.5*ruyy-1.5*rusq);
		eq[2][0] = w*(-3.*rux+4.5*ruxx-1.5*rusq);
		eq[3][0] = w*(-3.*ruy+4.5*ruyy-1.5*rusq);
	}

	/* v4 - v7 */
	{
		const double w = 1./36.;
		eq[4][0] = w*(3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		eq[5][0] = w*(3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
		eq[6][0] = w*(3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		eq[7][0] = w*(3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
	}
	
	/* v8 */
	{
		const double w = 4./9.;
		eq[8][0] = r + w*(-1.5*rusq);
	}
	
	/* phi eq dist. */
	/* v0 - v3 */
	{
		const double w = 1./9.;
		eq[0][1] = w*( 3.*pux+4.5*puxx-1.5*pusq);
		eq[1][1] = w*( 3.*puy+4.5*puyy-1.5*pusq);
		eq[2][1] = w*(-3.*pux+4.5*puxx-1.5*pusq);
		eq[3][1] = w*(-3.*puy+4.5*puyy-1.5*pusq);
	}

	/* v4 - v7 */
	{
		const double w = 1./36.;
		eq[4][1] = w*(3.*( pux+puy)+4.5*(puxx+puyy+puxy)-1.5*pusq);
		eq[5][1] = w*(3.*(-pux+puy)+4.5*(puxx+puyy-puxy)-1.5*pusq);
		eq[6][1] = w*(3.*(-pux-puy)+4.5*(puxx+puyy+puxy)-1.5*pusq);
		eq[7][1] = w*(3.*( pux-puy)+4.5*(puxx+puyy-puxy)-1.5*pusq);
	}
	
	/* v8 */
	{
		const double w = 4./9.;
		eq[8][1] = p + w*(-1.5*pusq);
	}
}

double 
lb_Landau_pressure_3d(const LBLandauParameters* params,
		     const double* LB_RESTRICT rho,
			 const double* LB_RESTRICT phi,
		     const double* LB_RESTRICT drx,
		     const double* LB_RESTRICT dry,
		     const double* LB_RESTRICT drz,
			 const double* LB_RESTRICT dpx,
		     const double* LB_RESTRICT dpy,
		     const double* LB_RESTRICT dpz,
			 const double* LB_RESTRICT ddr,
		     const double* LB_RESTRICT ddp)
{
	lb_assert(params && rho && phi && drx && dry && drz && dpx && dpy && dpz && ddr && ddp);

	return (1./3.)*(*rho) + 0.5*params->a*(*phi)*(*phi) + 0.75*params->b*(*phi)*(*phi)*(*phi)*(*phi) 
			- params->K*( (*phi)*(*ddp) + 0.5*( (*dpx)*(*dpx) + (*dpy)*(*dpy) + (*dpy)*(*dpy) + (*dpz)*(*dpz) ) );
}

void
lb_d3q19_Landau_equilibrium(const LBLandauParameters* params,
				const double* LB_RESTRICT rhop,
				const double* LB_RESTRICT up,
				const double* LB_RESTRICT drx,
				const double* LB_RESTRICT dry,
				const double* LB_RESTRICT drz,
				const double* LB_RESTRICT dpx,
				const double* LB_RESTRICT dpy,
				const double* LB_RESTRICT dpz,
				const double* LB_RESTRICT ddr,
				const double* LB_RESTRICT ddp,
				double** LB_RESTRICT eq)
{
	const double r = *rhop++;
	const double p = *rhop;

	const double ux = *up++;
	const double uy = *up++;
	const double uz = *up;
	
	const double mu = lb_Landau_chemical_potential(params,&r,&p,ddp);
	const double Pr = lb_Landau_pressure_3d(params,&r,&p,drx,dry,drz,dpx,dpy,dpz,ddr,ddp);
	
	const double muu = r*mu;
	const double Pxx = Pr + params->K*(*dpx)*(*dpx);
	const double Pyy = Pr + params->K*(*dpy)*(*dpy);
	const double Pzz = Pr + params->K*(*dpz)*(*dpz);
	const double Pxy =      params->K*(*dpx)*(*dpy);
	const double Pxz =      params->K*(*dpx)*(*dpz);
	const double Pyz =      params->K*(*dpy)*(*dpz);

	/* rho */
	const double rux = r*ux;
	const double ruy = r*uy;
	const double ruz = r*uz;
	const double ruxx = rux*ux + Pxx;
	const double ruyy = ruy*uy + Pyy;
	const double ruzz = ruz*uz + Pzz;	
	const double ruxy = 2.*rux*uy + 2.*Pxy;
	const double ruxz = 2.*rux*uz + 2.*Pxz;
	const double ruyz = 2.*ruy*uz + 2.*Pyz;	
	const double rusq = ruxx + ruyy + ruzz;

	/* phi */
	const double pux = p*ux;
	const double puy = p*uy;
	const double puz = p*uz;
	const double puxx = pux*ux + muu;
	const double puyy = puy*uy + muu;
	const double puzz = puz*uz + muu;	
	const double puxy = 2.*pux*uy;
	const double puxz = 2.*pux*uz;
	const double puyz = 2.*puy*uz;	
	const double pusq = puxx + puyy + puzz;
	
	/* rho eq dist */
	/* v0 to v5 */
	{
		const double w = 1./18.;
	
		eq[0][0] = w*( 3.*rux+4.5*ruxx-1.5*rusq);
		eq[1][0] = w*(-3.*rux+4.5*ruxx-1.5*rusq);

		eq[2][0] = w*( 3.*ruy+4.5*ruyy-1.5*rusq);
		eq[3][0] = w*(-3.*ruy+4.5*ruyy-1.5*rusq);

		eq[4][0] = w*( 3.*ruz+4.5*ruzz-1.5*rusq);
		eq[5][0] = w*(-3.*ruz+4.5*ruzz-1.5*rusq);
	}

	/* v6 to v17 */
	{
		const double w = 1./36.;
			
		/* z=0 plane */
		{
			eq[6][0] = w*(3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
			eq[7][0] = w*(3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
			
			eq[8][0] = w*(3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
			eq[9][0] = w*(3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		}

		/* y=0 plane */
		{
			eq[10][0] = w*(3.*( rux+ruz)+4.5*(ruxx+ruzz+ruxz)-1.5*rusq);
			eq[11][0] = w*(3.*(-rux+ruz)+4.5*(ruxx+ruzz-ruxz)-1.5*rusq);
			
			eq[12][0] = w*(3.*( rux-ruz)+4.5*(ruxx+ruzz-ruxz)-1.5*rusq);
			eq[13][0] = w*(3.*(-rux-ruz)+4.5*(ruxx+ruzz+ruxz)-1.5*rusq);
		}

		/* x=0 plane */
		{
			eq[14][0] = w*(3.*( ruy+ruz)+4.5*(ruyy+ruzz+ruyz)-1.5*rusq);
			eq[15][0] = w*(3.*(-ruy+ruz)+4.5*(ruyy+ruzz-ruyz)-1.5*rusq);
			
			eq[16][0] = w*(3.*( ruy-ruz)+4.5*(ruyy+ruzz-ruyz)-1.5*rusq);
			eq[17][0] = w*(3.*(-ruy-ruz)+4.5*(ruyy+ruzz+ruyz)-1.5*rusq);
		}
	}

	{
		const double w = 1./3.;
		eq[18][0] = r + w*(-1.5*rusq);
	}
	/* phi eq dist */
	/* v0 to v5 */
	{
		const double w = 1./18.;
	
		eq[0][1] = w*( 3.*pux+4.5*puxx-1.5*pusq);
		eq[1][1] = w*(-3.*pux+4.5*puxx-1.5*pusq);

		eq[2][1] = w*( 3.*puy+4.5*puyy-1.5*pusq);
		eq[3][1] = w*(-3.*puy+4.5*puyy-1.5*pusq);

		eq[4][1] = w*( 3.*puz+4.5*puzz-1.5*pusq);
		eq[5][1] = w*(-3.*puz+4.5*puzz-1.5*pusq);
	}

	/* v6 to v17 */
	{
		const double w = 1./36.;
			
		/* z=0 plane */
		{
			eq[6][1] = w*(3.*( pux+puy)+4.5*(puxx+puyy+puxy)-1.5*pusq);
			eq[7][1] = w*(3.*(-pux+puy)+4.5*(puxx+puyy-puxy)-1.5*pusq);
			
			eq[8][1] = w*(3.*( pux-puy)+4.5*(puxx+puyy-puxy)-1.5*pusq);
			eq[9][1] = w*(3.*(-pux-puy)+4.5*(puxx+puyy+puxy)-1.5*pusq);
		}

		/* y=0 plane */
		{
			eq[10][1] = w*(3.*( pux+puz)+4.5*(puxx+puzz+puxz)-1.5*pusq);
			eq[11][1] = w*(3.*(-pux+puz)+4.5*(puxx+puzz-puxz)-1.5*pusq);
			
			eq[12][1] = w*(3.*( pux-puz)+4.5*(puxx+puzz-puxz)-1.5*pusq);
			eq[13][1] = w*(3.*(-pux-puz)+4.5*(puxx+puzz+puxz)-1.5*pusq);
		}

		/* x=0 plane */
		{
			eq[14][1] = w*(3.*( puy+puz)+4.5*(puyy+puzz+puyz)-1.5*pusq);
			eq[15][1] = w*(3.*(-puy+puz)+4.5*(puyy+puzz-puyz)-1.5*pusq);
			
			eq[16][1] = w*(3.*( puy-puz)+4.5*(puyy+puzz-puyz)-1.5*pusq);
			eq[17][1] = w*(3.*(-puy-puz)+4.5*(puyy+puzz+puyz)-1.5*pusq);
		}
	}

	{
		const double w = 1./3.;
		eq[18][1] = p + w*(-1.5*pusq);
	}
}
