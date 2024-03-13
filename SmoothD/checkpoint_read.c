/***************************************************************

 checkpoint_read.c

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


/***************************************************************

  reads parameter space from binary file "checkpointXXXX.out"

***************************************************************/
#include "TP06_OpSplit_2D.h"

int checkpoint_read( double **u, double *time, int *t, int count, int N)
{
  int elements_to_read, i;
  int n, m, M;

  double u_n_m;
  char fname[80];

  FILE *chkpt_file;

  M = NUM_STATES;
  sprintf( fname,"%s%06d.out",CHKPTROOT,count );
  printf("opening checkpoint file %s\n",fname);
  chkpt_file = fopen( fname, "rb" );

  fread( time, sizeof(double),1,chkpt_file);
  fread( t, sizeof(int),1,chkpt_file);

  printf("read time = %g, t = %d\n",*time,*t);
  elements_to_read = N * M;
  printf("reading %d elements\n",elements_to_read);

  i = 0;
  u_n_m = 0.0;
  for (n = 1; n <= N; n++)
  {
	for (m = 1; m <= M; m++)
    {
      i += fread( &u_n_m, sizeof(double), 1, chkpt_file );
      u[n][m] = u_n_m;
    }
  }

  printf("read %d elements\n", i);
  if (ferror(chkpt_file)) perror("error reading data");

  fclose(chkpt_file);
  return (1);
}
