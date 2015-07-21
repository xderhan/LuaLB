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

#include "lb-Mix.h"
#include "lb-messages.h"

double lb_Mix_chemical_potential(const LBMixParameters* params,
			    const double* LB_RESTRICT rho,
			    const double* LB_RESTRICT phi,
			    const double* LB_RESTRICT ddp)
{
	lb_assert(params && rho && phi && ddp);

	return 0.5*params->T*log( (*rho + *phi)/(*rho - *phi) ) - params->K*(*ddp)
			- 0.5*params->lamda*(*phi)/(*rho);
}

double
lb_Mix_pressure_2d(const LBMixParameters* params,
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

	return params->T*(*rho)	-     params->K*( (*rho)*(*ddr) + (*phi)*(*ddp) ) 
							- 0.5*params->K*( (*drx)*(*drx) + (*dry)*(*dry) 
											+ (*dpx)*(*dpx) + (*dpy)*(*dpy) );
}

void
lb_d2q9_Mix_equilibrium(const LBMixParameters* params,
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
	const double rho = *averages++;
	const double phi = *averages;

	const double ux = *u++;
	const double uy = *u;

	const double rux = rho*ux;
	const double ruy = rho*uy;

	const double rux2 = rux*ux;
	const double ruy2 = ruy*uy;

	/* Phi */
	const double pux = phi*ux;
	const double puy = phi*uy;

	const double pux2 = pux*ux;
	const double puy2 = puy*uy;

	const double drx2 = params->K*(*drx)*(*drx);
	const double dry2 = params->K*(*dry)*(*dry);

	const double dpx2 = params->K*(*dpx)*(*dpx);
	const double dpy2 = params->K*(*dpy)*(*dpy);

	const double rddr = params->K*( rho*(*ddr) + phi*(*ddp) );

	/*const double brho = params->b*rho;
	const double brho1 = 1.0 - brho;
	const double Td = params->T/brho1;
	const double arho = params->a*rho;*/

	//const double p0 = rho*(Td - arho);
	//const double trP = (p0 - rddr)/3;
	//const double trP = (p0 - rddr)/8.;

	//const double dp0 = Td*(1.0 + brho/brho1) - 2*arho;
	//const double lamda = (params->tau0 - 0.5)*(1.0/3.0 - dp0);
	//const double nuudr = nu*(ux*(*drx) + uy*(*dry));

	const double trP = (drx2 + dry2 + dpx2 + dpy2)/16. + lb_Mix_pressure_2d(params, 
														&rho, &phi,
														drx, dry, 
														dpx, dpy,
														ddr, ddp)/8.;

	const double nu = lb_Mix_chemical_potential(params, &rho, &phi, ddp);

	const double Gxx = params->K*( (*dpx)*(*dpx) - (*dpy)*(*dpy) )/16.;
	const double Gyy = -Gxx;

	const double Gxy = params->K*(*dpx)*(*dpy)/8.;
	const double Gyx = Gxy;	

	const double H = params->Gr*nu/8.;

	{
		const double t1 = 2.*trP;
		const double t2 = t1 - (rux2 + ruy2)/8.;

		const double t3 = rux/3.;
		const double t4 = ruy/3.;

		const double t5 = 0.5*rux2;
		const double t6 = 0.5*ruy2;

		const double t7 = 4.*Gxx;
		const double t8 = 4.*Gyy;

		eq[0][0] = t2 + t3 + t5 + t7;
		eq[1][0] = t2 + t4 + t6 + t8;
		eq[2][0] = t2 - t3 + t5 + t7;
		eq[3][0] = t2 - t4 + t6 + t8;

	}

	{
		const double t1 = trP;
		const double t2 = t1 + (rux2 + ruy2)/16.;

		const double t3 = rux/12.;
		const double t4 = ruy/12.;

		const double t5 = 0.25*rux*uy;

		const double t6 = Gxy;
		const double t7 = Gyx;

		const double t8 = t5 + t6 + t7;

		eq[4][0] = t2 + t3 + t4 + t8;
		eq[5][0] = t2 - t3 + t4 - t8;
		eq[6][0] = t2 - t3 - t4 + t8;
		eq[7][0] = t2 + t3 - t4 - t8;

	}

	eq[8][0] = rho - 12.*trP - 3.*(rux2 + ruy2)/4.;

	/* Now for g_eq */
	{
		const double t1 = 2.*H;
		const double t2 = t1 - (pux2 + puy2)/8.;

		const double t3 = pux/3.;
		const double t4 = puy/3.;

		const double t5 = 0.5*pux2;
		const double t6 = 0.5*puy2;

		eq[0][1] = t2 + t3 + t5;
		eq[1][1] = t2 + t4 + t6;
		eq[2][1] = t2 - t3 + t5;
		eq[3][1] = t2 - t4 + t6;
	}

	{
		const double t1 = H;
		const double t2 = t1 + (pux2 + puy2)/16.;

		const double t3 = pux/12.;
		const double t4 = puy/12.;

		const double t5 = 0.25*pux*uy;

		eq[4][1] = t2 + t3 + t4 + t5;
		eq[5][1] = t2 - t3 + t4 - t5;
		eq[6][1] = t2 - t3 - t4 + t5;
		eq[7][1] = t2 + t3 - t4 - t5;
	}

	eq[8][1] = phi - 12.*H - 3.*(pux2 + puy2)/4.;

	/* for testing only
		eq[0][0] = 0.1;
		eq[1][0] = 0.;
		eq[2][0] = 0.1;
		eq[3][0] = 0.;

		eq[4][0] = 0.1;
		eq[5][0] = 0.;
		eq[6][0] = 0.;
		eq[7][0] = 0.;

		eq[0][1] = 0.;
		eq[1][1] = 0.1;
		eq[2][1] = 0.;
		eq[3][1] = 0.;

		eq[4][1] = 0.;
		eq[5][1] = 0.;
		eq[6][1] = 0.;
		eq[7][1] = 0.1; */
}

double 
lb_Mix_pressure_3d(const LBMixParameters* params,
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

	return params->T*(*rho)	-     params->K*( (*rho)*(*ddr) + (*phi)*(*ddp) ) 
							- 0.5*params->K*( (*drx)*(*drx) + (*dry)*(*dry) + (*drz)*(*drz)
											+ (*dpx)*(*dpx) + (*dpy)*(*dpy) + (*dpz)*(*dpz) );
}

void
lb_d3q19_Mix_equilibrium(const LBMixParameters* params,
				const double* LB_RESTRICT rhop,
				const double* LB_RESTRICT up,
				const double* LB_RESTRICT drx,
				const double* LB_RESTRICT dry,
				const double* LB_RESTRICT drz,
				const double* LB_RESTRICT dpx,
				const double* LB_RESTRICT dpy,
				const double* LB_RESTRICT dpz,
				const double* LB_RESTRICT ddr,
				const double* LB_RESTRICT ddrp,
				double** LB_RESTRICT eq)
{
	const double rho = *rhop++;
	const double phi = *rhop;

	const double ux = *up++;
	const double uy = *up++;
	const double uz = *up;

	const double rux = rho*ux;
	const double ruy = rho*uy;
	const double ruz = rho*uz;

	const double rux2 = rux*ux;
	const double ruy2 = ruy*uy;
	const double ruz2 = ruz*uz;

	/* Phi */
	const double pux = phi*ux;
	const double puy = phi*uy;
	const double puz = phi*uz;

	const double pux2 = pux*ux;
	const double puy2 = puy*uy;
	const double puz2 = puz*uz;

	const double drx2 = params->K*(*drx)*(*drx);
	const double dry2 = params->K*(*dry)*(*dry);
	const double drz2 = params->K*(*drz)*(*drz);

	const double dpx2 = params->K*(*dpx)*(*dpx);
	const double dpy2 = params->K*(*dpy)*(*dpy);
	const double dpz2 = params->K*(*dpz)*(*dpz);

	/*const double brho = params->b*rho;
	const double brho1 = 1.0 - brho;
	const double Td = params->T/brho1;
	const double arho = params->a*rho;
	const double p0 = rho*(Td - arho);

	const double Pcc = p0 - params->K*rho*(*ddr);
	const double Pxx = Pcc + 0.5*(drx2 - dry2 - drz2);
	const double Pxy = params->K*(*drx)*(*dry);
	const double Pxz = params->K*(*drx)*(*drz);
	const double Pyy = Pcc + 0.5*(dry2 - drx2 - drz2);
	const double Pyz = params->K*(*dry)*(*drz);
	const double Pzz = Pcc + 0.5*(drz2 - drx2 - dry2);

	const double trP = (Pxx + Pyy + Pzz)/9.0;
	const double ru2 = rux2 + ruy2 + ruz2;

	const double dp0 = Td*(1.0 + brho/brho1) - 2*arho;
	const double nu = (params->tau0 - 0.5)*(1.0/3.0 - dp0);

	const double ndrx = nu*(*drx);
	const double ndry = nu*(*dry);
	const double ndrz = nu*(*drz);

	const double nuudr = ux*ndrx + uy*ndry + uz*ndrz;

	{
		const double t1 = 0.25*(ru2 - trP) - nuudr;

		{
			const double t2 = t1 + 0.25*(Pxx + rux2) + 0.5*ndrx*ux;
			const double t3 = rux/6.0;

			eq[0][0] = t2 + t3;
			eq[1][0] = t2 - t3;
		}

		{
			const double t2 = t1 + 0.25*(Pyy + ruy2) + 0.5*ndry*uy;
			const double t3 = ruy/6.0;

			eq[2][0] = t2 + t3;
			eq[3][0] = t2 - t3;
		}

		{
			const double t2 = t1 + 0.25*(Pzz + ruz2) + 0.5*ndrz*uz;
			const double t3 = ruz/6.0;

			eq[4][0] = t2 + t3;
			eq[5][0] = t2 - t3;
		}
	}

	{
		const double t1 = 0.5*(0.5*nuudr - trP - 0.25*ru2);

		{
			const double t2 = t1 + 0.125*(Pxx + Pyy);

			{
				const double t3 = 0.25*(ux + uy);
				const double t4 = rux + ruy;
				const double t5 = t2 + 0.5*t3*t4 
						+ t3*(ndrx + ndry) + 0.25*Pxy;
				const double t6 = t4/12.0;

				eq[6][0] = t5 + t6;
				eq[9][0] = t5 - t6;
			}

			{
				const double t3 = 0.25*(ux - uy);
				const double t4 = rux - ruy;
				const double t5 = t2 + 0.5*t3*t4
						+ t3*(ndrx - ndry) - 0.25*Pxy;
				const double t6 = t4/12.0;

				eq[7][0] = t5 - t6;
				eq[8][0] = t5 + t6;
			}
		}

		{
			const double t2 = t1 + 0.125*(Pxx + Pzz);

			{
				const double t3 = 0.25*(ux + uz);
				const double t4 = rux + ruz;
				const double t5 = t2 + 0.5*t3*t4 
						+ t3*(ndrx + ndrz) + 0.25*Pxz;
				const double t6 = t4/12.0;

				eq[10][0] = t5 + t6;
				eq[13][0] = t5 - t6;
			}

			{
				const double t3 = 0.25*(ux - uz);
				const double t4 = rux - ruz;
				const double t5 = t2 + 0.5*t3*t4
						+ t3*(ndrx - ndrz) - 0.25*Pxz;
				const double t6 = t4/12.0;

				eq[11][0] = t5 - t6;
				eq[12][0] = t5 + t6;
			}
		}

		{
			const double t2 = t1 + 0.125*(Pyy + Pzz);

			{
				const double t3 = 0.25*(uy + uz);
				const double t4 = ruy + ruz;
				const double t5 = t2 + 0.5*t3*t4 
						+ t3*(ndry + ndrz) + 0.25*Pyz;
				const double t6 = t4/12.0;

				eq[14][0] = t5 + t6;
				eq[17][0] = t5 - t6;
			}

			{
				const double t3 = 0.25*(uy - uz);
				const double t4 = ruy - ruz;
				const double t5 = t2 + 0.5*t3*t4
						+ t3*(ndry - ndrz) - 0.25*Pyz;
				const double t6 = t4/12.0;

				eq[15][0] = t5 - t6;
				eq[16][0] = t5 + t6;
			}
		}
	}

	eq[18][0] = rho - 6.0*trP - 1.5*ru2;*/
}
