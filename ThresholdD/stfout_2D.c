/***************************************************************

 stfout.c

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

/***************************************************************

  stfout

  prints full parameter space u to an stf file

***************************************************************/

int stfout_2D( double **u, int **geomarray, int stfcount, int nx, int ny )
{
  int lay, row, col;
  int index, outint;
  double outdouble;
  char fname[80];
  char gzipfile[80];

  FILE *stf_file;

  sprintf( fname,"%s%04d.stf",STFFILEROOT,stfcount );

  stf_file = fopen( fname, "w" );

  fprintf(stf_file, "NAME Vm\n");
  fprintf(stf_file, "RANK 2\n");
  fprintf(stf_file, "DIMENSIONS %d %d\n",nx, ny);
  fprintf(stf_file, "BOUNDS %d %d %d %d\n", 0, nx-1, 0, ny-1);
  fprintf(stf_file, "SCALAR\n");
  fprintf(stf_file, "DATA\n");

  for (row = 1;row <= ny; row+=1)
    {
    for (col = 1; col <= nx; col+=1)
      {
      index = (geomarray[row][col] > 0)?geomarray[row][col]:0;
      if (index > 0)
          outdouble = u[index][1];
      else
          outdouble = -100.0;
      fprintf(stf_file, "%4.2f ", outdouble);
      }
    fprintf( stf_file, "\n");
  }
  fclose( stf_file);

  sprintf( gzipfile,"gzip %s",fname );
  if (system(gzipfile) == -1 ) nrerror("failed to gzip file\n");

  return (1);
}
