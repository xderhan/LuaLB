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

#ifndef LB_D3Q19_DENSITY_GUARDS_H
#define LB_D3Q19_DENSITY_GUARDS_H

/*
 *  config.h should be already included in order to have proper LB_ENABLE_MPI
 */

#ifdef LB_ENABLE_MPI
#include <mpi.h>
#endif /* LB_ENABLE_MPI */

#include <lb-macros.h>
#include <lb-partition-info.h>

#include <lb-d3q19.h>

LB_BEGIN_DECLS



void

lb_d3q19_average_pdfs(LBD3Q19*, const double **, double *, double *);

void
lb_d3q19_setup_rho_guards(LBD3Q19*, double*);

LB_END_DECLS

#endif /* LB_D3Q19_DENSITY_GUARDS_H */
