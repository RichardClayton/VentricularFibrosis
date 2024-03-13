
/***************************************************************

 initialise_variables.c

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

#include <TP06_OpSplit_2D.h>

void initialise_variables_2D( double **u, int N )
{

  int n,m;

  /* Indices for u array */
  int V =      1;
  int M =      2;
  int H =      3;
  int J =      4;
  int R =      5;
  int S =      6;
  int D =      7;
  int F =      8;
  int F2 =     9;
  int FCass = 10;
  int Xr1 =   11;
  int Xr2 =   12;
  int Xs =    13;
  int RR =    14;
  int OO =    15;
  int CaSS =  16;
  int CaSR =  17;
  int Cai =   18;
  int Nai =   19;
  int Ki =    20;
  
  /* Initial Gate Conditions */
  // original commented out
  // new values from CellML
  for (n = 1; n <= N; n++)
    {
  /*  u[n][V] = -86.2;
    u[n][M] = 0.0;
    u[n][H] = 0.75;
    u[n][J] = 0.75;
    u[n][Xr1] = 0.0;
    u[n][Xr2] = 1.0;
    u[n][Xs] = 0.0;
    u[n][R] = 0.0;
    u[n][S] = 1.0;
    u[n][D] = 0.0;
    u[n][F] = 1.0;
    u[n][F2] = 1.0;
    u[n][FCass] = 1.0;
    u[n][RR] = 1.0;
    u[n][OO] = 0.0;
    u[n][Cai] = 0.00007;
    u[n][CaSR] = 3.0; //1.3;
    u[n][CaSS] = 0.00007;
    u[n][Nai] = 7.67;
    u[n][Ki] = 138.3; */

    u[n][V]     = -85.23;
    u[n][M]     = 0.00172;
    u[n][H]     = 0.7444;
    u[n][J]     = 0.7045;
    u[n][Xr1]   = 0.000621;
    u[n][Xr2]   = 0.4712;
    u[n][Xs]    = 0.0095;
    u[n][R]     = 0.0000000242;
    u[n][S]     = 0.999998;
    u[n][D]    = 0.00003373;
    u[n][F]     = 0.7888;
    u[n][F2]    = 0.9755;
    u[n][FCass] = 0.9953;
    u[n][RR]    = 0.9073;
    u[n][OO]    = 0.0;
    u[n][Cai]   = 0.000126;
    u[n][CaSR]  = 3.64;
    u[n][CaSS]  = 0.00036;
    u[n][Nai]   = 8.604;
    u[n][Ki]    = 136.89;

    }

}
