#ifndef LB_BGK_H
#define LB_BGK_H

#include "lb-macros.h"

LB_BEGIN_DECLS

typedef struct {
	double T;
	double a;
	double b;
	double K;
	double tau;
} LBGKParameters;

double lb_BGK_pressure_2d(const LBGKParameters*,
			    const double* rho,
			    const double* drx,
			    const double* dry,
			    const double* ddr);

void lb_d2q9_BGK_equilibrium(const LBGKParameters*,
			       const double* rho,
				   const double* u,
			       double equilibrium_pdfs[9]);

double lb_BGK_pressure_3d(const LBGKParameters*,
			    const double* rho,
			    const double* drx,
			    const double* dry,
			    const double* drz,
			    const double* ddr);

void lb_d3q19_BGK_equilibrium(const LBGKParameters*,
			        const double* rho,
					const double* u,
			        double equilibrium_pdfs[19]);

LB_END_DECLS

#endif /* LB_BGK_H */
