
#include "TP06_OpSplit_2D.h"

double diffusion_2D_modD(double **u, int **nneighb, int n, int N, double *D, double dx2)
{
/* Work out isotropic diffusion */

  int V = 1;

  double diffusion;
  double d2vdx2, d2vdy2;
  double nn6, nn2, nn8, nn4;
  double Dnn2, Dnn4, Dnn6, Dnn8;
  double dvdx, dvdy,dDdx,dDdy;
  double twodx;

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

  twodx = DX * 2.0;

  dvdx = (nn6 - nn2) / twodx;
  dvdy = (nn8 - nn4) / twodx;

  dDdx = ((Dnn6 > 0) && (Dnn2 > 0) && (D[n] > 0)) ? (Dnn6 - Dnn2) / twodx : 0.0;
  dDdy = ((Dnn8 > 0) && (Dnn4 > 0) && (D[n] > 0)) ? (Dnn8 - Dnn4) / twodx : 0.0;

  d2vdx2 = (nn6 + nn2 - (2.0 * u[n][V]))/dx2;
  d2vdy2 = (nn4 + nn8 - (2.0 * u[n][V]))/dx2;

  //d2vdx2 = (((Dnn6 + D[n])*(nn6 - u[n][V])) - ((Dnn2 + D[n])*(u[n][V] - nn2)))/(2.0 * dx2);
  //d2vdy2 = (((Dnn4 + D[n])*(nn4 - u[n][V])) - ((Dnn8 + D[n])*(u[n][V] - nn8)))/(2.0 * dx2);

  diffusion = D[n] * (d2vdx2 + d2vdy2) + dvdx*dDdx + dvdy*dDdy;

  return(diffusion);

}
