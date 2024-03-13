
/***************************************************************

 diffusion_2D_modD.c

 Version       2.1

 Date          23-oct-2023

 Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)
 
 This file is part of VentricularFibrosis.

 Copyright (c) Richard Clayton, 
 Department of Computer Science, 
 University of Sheffield, 2006, 2009, 2016, 2023

 VentricularFibrosis is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*********************************************************************/

#include "TP06_OpSplit_2D.h"

double diffusion_2D_modD(double **u, int **nneighb, int n, int N, double *D, double dx2)
{
/* Work out isotropic diffusion with weighted diffusion coefficient */

  int V = 1;

  double diffusion;
  double d2vdx2, d2vdy2;
  double nn6, nn2, nn8, nn4;
  double Dnn2, Dnn4, Dnn6, Dnn8;

  /* Calculate Vm of nearest neighbours, allowing for no-flux
     at the boundaries */

  nn6 = (nneighb[n][6] > 0) ? u[nneighb[n][6]][V] : u[n][V];
  nn2 = (nneighb[n][2] > 0) ? u[nneighb[n][2]][V] : u[n][V];
  nn4 = (nneighb[n][4] > 0) ? u[nneighb[n][4]][V] : u[n][V];
  nn8 = (nneighb[n][8] > 0) ? u[nneighb[n][8]][V] : u[n][V];

  Dnn6 = (nneighb[n][6] > 0) ? D[nneighb[n][6]] : D[n];
  Dnn2 = (nneighb[n][2] > 0) ? D[nneighb[n][2]] : D[n];
  Dnn4 = (nneighb[n][4] > 0) ? D[nneighb[n][4]] : D[n];
  Dnn8 = (nneighb[n][8] > 0) ? D[nneighb[n][8]] : D[n];

  //d2vdx2 = (nn6 + nn2 - (2.0 * u[n][V]))/dx2;
  //d2vdy2 = (nn4 + nn8 - (2.0 * u[n][V]))/dx2;

  d2vdx2 = (((Dnn6 + D[n])*(nn6 - u[n][V])) - ((Dnn2 + D[n])*(u[n][V] - nn2)))/(2.0 * dx2);
  d2vdy2 = (((Dnn4 + D[n])*(nn4 - u[n][V])) - ((Dnn8 + D[n])*(u[n][V] - nn8)))/(2.0 * dx2);

  diffusion = d2vdx2 + d2vdy2;

  return(diffusion);

}
