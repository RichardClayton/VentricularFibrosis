/***************************************************************

 initialise_geometry.c

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

/* this version sets up smoothly varying diffusion with no-flux boundaries */

int initialise_geometry_2D( int **geom, int nrows, int ncols, int **nneighb, double *D )
{

  int row, col, n, m;
  int N;
  float dTemp = 0.0;
  FILE *inFile;

  /* read in diffusion information from file */
  n = 1;
  inFile = fopen(DIFFUSIONFILE,"r");
  for (row = 1; row <= nrows; row++)
    {
    for (col = 1; col <= ncols; col++)
      {
      if (fscanf(inFile, "%f ", &dTemp) == 0)
        nrerror("error reading DIFFUSIONFILE\n");

      if (dTemp <= 0.025)
        {
        geom[row][col] = -1;
        }

      else
        {
        geom[row][col] = n;
        D[n] = (double) dTemp;
        n++;
        }
      }
    }

  printf("Successfully read %d entries from DIFFUSIONFILE\n",n);
  fclose(inFile);
  N = n;

  /* now work out nearest neighbours */

  /*  1 2 3
      8   4
      7 6 5  */

  /* grid points on the boundary are identified by nneighb of 0 */

  for (row = 1; row <= nrows; row++ ) for (col = 1; col <= ncols; col++)
    {
    if (geom[row][col] > 0)
      {
      n = geom[row][col];
      nneighb[n][1] = ((row>1)&&(col>1))         ? geom[row-1][col-1] :0;
      nneighb[n][2] = (row>1)                    ? geom[row-1][col]   :0;
      nneighb[n][3] = ((row>1)&&(col<ncols))     ? geom[row-1][col+1] :0;
      nneighb[n][4] = (col<ncols)                ? geom[row][col+1]   :0;
      nneighb[n][5] = ((row<nrows)&&(col<ncols)) ? geom[row+1][col+1] :0;
      nneighb[n][6] = (row<nrows)                ? geom[row+1][col]   :0;
      nneighb[n][7] = ((row<nrows)&&(col>1))     ? geom[row+1][col-1] :0;
      nneighb[n][8] = (col>1)                    ? geom[row][col-1]   :0;

      /* test allocation */
      //printf("r %d, c %d, n%d -- %d %d %d %d\n",row,col,n,nneighb[n][2],nneighb[n][4],nneighb[n][6],nneighb[n][8]);

      /* implement ring geometry */
      //nneighb[n][4] = (col<ncols)                ? geom[row][col+1]   :geom[row][1]; //0;
      //nneighb[n][8] = (col>1)                    ? geom[row][col-1]   :geom[row][ncols];//0;
      //n++;
      }
    }
    return(N);
}
