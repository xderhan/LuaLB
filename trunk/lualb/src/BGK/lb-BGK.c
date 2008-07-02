#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "lb-BGK.h"
#include "lb-messages.h"

double
lb_BGK_pressure_2d(const LBGKParameters* params,
		     const double* LB_RESTRICT rho,
		     const double* LB_RESTRICT drx,
		     const double* LB_RESTRICT dry,
		     const double* LB_RESTRICT ddr)
{
	lb_assert(params && rho && drx && dry && ddr);

	return (*rho)*(params->T/(1.0 - params->b*(*rho))
		     - params->a*(*rho))
		- params->K*((*drx)*(*drx) + (*dry)*(*dry) + 0.5*(*rho)*(*ddr));
}

void
lb_d2q9_BGK_equilibrium(const LBGKParameters* params,
			const double* LB_RESTRICT rho,
			const double* LB_RESTRICT u,
			double* LB_RESTRICT eq)
{
	const double r = *rho;
	
	const double ux = *u++;
	const double uy = *u;

	const double rux = r*ux;
	const double ruy = r*uy;
	
	const double ruxx = rux*ux;
	const double ruyy = ruy*uy;
	
	const double ruxy = 2.*rux*uy;
	
	const double rusq = ruxx+ruyy;
	
	/* v0 - v3 */
	{
		const double w = 1./9.;
		*eq++ = w*(r+3.*rux+4.5*ruxx-1.5*rusq);
		*eq++ = w*(r+3.*ruy+4.5*ruyy-1.5*rusq);
		*eq++ = w*(r-3.*rux+4.5*ruxx-1.5*rusq);
		*eq++ = w*(r-3.*ruy+4.5*ruyy-1.5*rusq);
	}

	/* v4 - v7 */
	{
		const double w = 1./36.;
		*eq++ = w*(r+3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		*eq++ = w*(r+3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
		*eq++ = w*(r+3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		*eq++ = w*(r+3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
	}
	
	/* v8 */
	{
		const double w = 4./9.;
		*eq = w*(r-1.5*rusq);
	}
}

double
lb_BGK_pressure_3d(const LBGKParameters* params,
		     const double* LB_RESTRICT rho,
		     const double* LB_RESTRICT drx,
		     const double* LB_RESTRICT dry,
		     const double* LB_RESTRICT drz,
		     const double* LB_RESTRICT ddr)
{
	lb_assert(params && rho && drx && dry && drz && ddr);

	return (*rho)*(params->T/(1.0 - params->b*(*rho))
		     - params->a*(*rho))
		- params->K*((*drx)*(*drx) + (*dry)*(*dry) + (*drz)*(*drz)
			    + 0.5*(*rho)*(*ddr));
}

void
lb_d3q19_BGK_equilibrium(const LBGKParameters* params,
			   const double* LB_RESTRICT rhop,
			   const double* LB_RESTRICT up,
			   double* LB_RESTRICT eq)
{
	const double r = *rhop;

	const double ux = *up++;
	const double uy = *up++;
	const double uz = *up;

	const double rux = r*ux;
	const double ruy = r*uy;
	const double ruz = r*uz;

	const double ruxx = rux*ux;
	const double ruyy = ruy*uy;
	const double ruzz = ruz*uz;
	
	const double ruxy = 2.*rux*uy;
	const double ruxz = 2.*rux*uz;
	const double ruyz = 2.*ruy*uz;
	
	const double rusq = ruxx+ruyy+ruzz;

	/* v0 to v5 */
	{
		const double w = 1./18.;
	
		eq[0] = w*(r+3.*rux+4.5*ruxx-1.5*rusq);
		eq[1] = w*(r-3.*rux+4.5*ruxx-1.5*rusq);

		eq[2] = w*(r+3.*ruy+4.5*ruyy-1.5*rusq);
		eq[3] = w*(r-3.*ruy+4.5*ruyy-1.5*rusq);

		eq[4] = w*(r+3.*ruz+4.5*ruzz-1.5*rusq);
		eq[5] = w*(r-3.*ruz+4.5*ruzz-1.5*rusq);
	}

	/* v6 to v17 */
	{
		const double w = 1./36.;
			
		/* z=0 plane */
		{
			eq[6] = w*(r+3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
			eq[7] = w*(r+3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
			
			eq[8] = w*(r+3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq);
			eq[9] = w*(r+3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq);
		}

		/* y=0 plane */
		{
			eq[10] = w*(r+3.*( rux+ruz)+4.5*(ruxx+ruzz+ruxz)-1.5*rusq);
			eq[11] = w*(r+3.*(-rux+ruz)+4.5*(ruxx+ruzz-ruxz)-1.5*rusq);
			
			eq[12] = w*(r+3.*( rux-ruz)+4.5*(ruxx+ruzz-ruxz)-1.5*rusq);
			eq[13] = w*(r+3.*(-rux-ruz)+4.5*(ruxx+ruzz+ruxz)-1.5*rusq);
		}

		/* x=0 plane */
		{
			eq[14] = w*(r+3.*( ruy+ruz)+4.5*(ruyy+ruzz+ruyz)-1.5*rusq);
			eq[15] = w*(r+3.*(-ruy+ruz)+4.5*(ruyy+ruzz-ruyz)-1.5*rusq);
			
			eq[16] = w*(r+3.*( ruy-ruz)+4.5*(ruyy+ruzz-ruyz)-1.5*rusq);
			eq[17] = w*(r+3.*(-ruy-ruz)+4.5*(ruyy+ruzz+ruyz)-1.5*rusq);
		}
	}

	{
		const double w = 1./3.;
		eq[18] = w*(r-1.5*rusq);
	}
}
