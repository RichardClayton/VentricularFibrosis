/********************************************************************

 TP06_OpSplit_2D_main.c

 Version       2.1

 Date          23-oct-2023

 Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)

 This program evaluates the TP06 model equations for a
 2D sheet, with rapid pacing.

 The parameters for numerical solution are contained in the header file
 The PDE is solved using operator splitting, with a maximum dt
 of 10 x DT as defined in the header file.
  
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

void writeData(char *fname, double *dataToWrite, int **geom, int nrows, int ncols);

int main(int argc, char **argv)
{
  /* constants */
  const int tmax = NUM_ITERATIONS;			    // duration of simulation
  const int num_states = NUM_STATES;        // number of states stored at each grid point
  const int num_params = NUM_PARAMS;        // number of model parameters (inputs)
  const int voltage_steps = VOLTAGE_STEPS;  // voltage steps in the lookup table
  const int num_lookup = NUM_LOOKUP;		    // number of variables in lookup table
  const int nrows = ROWS;
  const int ncols = COLUMNS;
  const int RC = ROWS * COLUMNS;            // total number of grid points
  const int V = 1;                          // index of membrane voltage

  const double dtlong = DT;					        // long time step for diffusion
  const double half_dtlong = DT/2.0;		    // half time step for diffusion
  //const double D = DIFFUSION;				      // diffusion coefficient
  const double dx2 = DX*DX;					        // dx squared

  /* variables */
  int N;
  int **geom, **nneighb;    				        // arrays to store geometry and nearest neighbours
  int t, n, m, dummy;						            // array indices
  int k, ko, kmax;					                // parameters for adaptive timestep
  int row, col;								              // array indices
  int *celltype;							              // specify myocyte (1) or fibroblast (0)
  int stfcount = 0;						              // index for stf output

  double dtshort;							              // adaptive short time step for ODE solution
  double **u, **lookup;						          // arrays for storing model state and lookup table
  double time = 0.0;
  double timems = 0.0;
  double stimCurrent = 0.0;
  double *dVdt, *new_Vm, *old_Vm, *U;		    // arrays for storing state during updates
  double dummy1, dummy2;
  double dV;
  double *params;
  double *D;                               // array to hold local diffusion coefficient

  const double bcl = S1BCL;                // basic cycle length for pacing
  const int numS1Beats = NUMS1BEATS;       // number of S1 stimuli
  const int numS2Beats = NUMS2BEATS;       // number of extra beats delivered after last S1
  //const double s1s2 = S1S2;                // coupling interval to s2 stimulus
  //const double decrement = DECREMENT;      // decrement in pacing interval with each stimulus
  int s1Beat = 1;
  int s2Beat = 1;
  int radius = 5;                          // radius of stimulating electrode
  int radius2 = radius*radius;
  double nextStim = 0.0;                   // time of next stimulus in ms
  int S1stimFlag = 0;
  int S2stimFlag = 0;
  int n_75_75 = 0;                         // node of stimulus point
  int *rowList;                            // store rows indexed by n
  int *colList;                            // store cols indexed by n
  
  const double threshold = -70.0;          // threshold for APD90 detection
  const double lastS1 = bcl * (numS1Beats - 1.0);
  double **upStrokeTime;                   // array to store upstroke times for final S1 beat
  double **downStrokeTime;                 // array to store upstroke times for S2 beat
  double *timing;
  int *beat;

  FILE *egPtr;
  char outputFile[80];					          // filename for outputs

  /* Create geometry and nearest neighbour arrays */
  geom = imatrix( 1, nrows, 1, ncols );
  nneighb = imatrix( 1, RC, 1, 8 );
  printf("initialising geometry arrays\n");
  D = fvector(1, RC);
  N = initialise_geometry_2D( geom, nrows, ncols, nneighb, D);

  /* Initialise arrays */
  u = fmatrix( 1, N, 1, num_states );
  lookup = fmatrix( 0, num_lookup, 0, voltage_steps );
  dVdt = fvector( 1, N );
  new_Vm = fvector( 1, N );
  old_Vm = fvector( 1, N );
  U = fvector(1, num_states);
  params = fvector(1, num_params);
  celltype = ivector(1, N);

  /* arrays for apd90 detection */
  upStrokeTime = fmatrix(1, N, 1, numS1Beats + numS2Beats);
  downStrokeTime = fmatrix(1, N, 1, numS1Beats + numS2Beats);
  beat = ivector(1, N);
  timing = fvector(1,N);

  /* open files for output */
  sprintf(outputFile,"%sVm_2.txt",OUTPUTFILEROOT);
  egPtr = fopen(outputFile,"w");

  /* Initialise model state and paramaters */
  printf("initialising ...\n");
  initialise_variables_2D( u, N );
  printf("done\n");

  /* set celltype to be 1 throughout */
  m = 0;
  for (n = 1; n <= N; n++)
      celltype[n] = 1;

  /* initialise indexing of rows and columns */
  rowList = ivector(1,N);
  colList = ivector(1,N);
  n_75_75 = geom[75][75];

  for (row = 1; row <= nrows; row++)
    {
      for (col = 1; col <= ncols; col++)
        {
          n = geom[row][col];
          if (n >= 1)
            {
              rowList[n] = row;
              colList[n] = col;
              printf("n %d, row %d, col %d\n", n, row, col);
            }
        }
    }

  /* initialise upstroke and downstroke arrays */
  for (n = 1; n <= N; n++)
    {
      beat[n] = 1;
      for (m = 1; m <= numS1Beats + numS2Beats; m++)
        {
          upStrokeTime[n][m] = -1.0;
          downStrokeTime[n][m] = -1.0;
        }
    }

  /* Create lookup table */
  printf("create lookup table ...\n");
  dummy = create_TP06_lookup_OpSplit_2D( lookup );
  printf("done\n");

  /* last bit of initialisation */
  t = 0;
  time = 0.0;

  if (CHKPT_READ)
    {
  	dummy = checkpoint_read( u, &time, &t, CHKPT_READ_TIME, N );
  	stfcount = ceil(time);
  	t = t + 1;
    }

/******************************************/
/*               Main loop                */
/******************************************/
  printf("entering main loop\n");

  while (t < tmax)
	  {
      time = t * dtlong;
      t++;

/* step 1 */
/* only at start */
      if (t == 1)
         {
         for (n = 1; n <= N; n++)
            {
            new_Vm[n] = u[n][V] + half_dtlong * diffusion_2D_modD( u, nneighb, n, N, D, dx2 );
            dVdt[n] = new_Vm[n] - u[n][V];
            }
         }

/* step 2 */
/* set up integration with adaptive timestep */
      S1stimFlag = 0;
      S2stimFlag = 0;
      for (n = 1; n <= N; n++)
        {
          /* pacing protocol -- deliver stimulus to one corner */
          stimCurrent = 0.0;

          // S1 pacing
          if ((n == 1) && ((time <= 2.0) || ((time > 400.0)&&(time <= 402.0)) || ((time > 800.0)&&(time <= 802.0)))) // || ((time > 1200.0)&&(time <= 1201.0))))
            {
              printf("Preparing to deliver S1 stimulus at time %f, u[%d][1] = %f\n",time,n_75_75,u[n_75_75][1]);
              S1stimFlag = 1;
            }

          if ((n == 1) && (time > 900) && (u[n_75_75][1] <= -84.5) && (time < 2100))
            {
              printf("Preparing to deliver S2 stimulus at time %f, u[%d][1] = %f\n",time,n_75_75,u[n_75_75][1]);
              nextStim = time;
            }

          col = colList[n];
          col -= 75;
          row = rowList[n];
          row -= 75;

          if ((S1stimFlag == 1) && (row*row + col*col <= radius2))
            {
              stimCurrent = -52.0;
              printf("Delivering S1 -- time = %f, row = %d, col = %d, n = %d\n",time,row,col,n);

            }

          if ((nextStim > 0) && (time < nextStim + 2.0) && (time >= nextStim) && (row*row + col*col <= radius2))
            {
              printf("Delivering S2 -- time = %f, nextStim = %f\n",time,nextStim);
              stimCurrent = -52.0;
            }

          if ((nextStim > 0) && (time >= nextStim + 2.0))
            {
              nextStim = 0;
            }      

          /* Operator splitting with adaptive time step for ODE */
      	  u[n][V] = new_Vm[n];

          // uncomment these lines to implement adaptive time step
          // this implementation provides good agreement with standard scheme for dt=0.01 ms
          // apart from delay of ~0.1 ms in onset of AP upstroke
   	      if (dVdt[n] > 0.01) ko = 5; else ko = 1;
              kmax = ko + floor(fabs(dVdt[n]) * 20.0); // kmax varies between ko and 10
          if (kmax > ceil(dtlong/0.01))
              kmax = dtlong/0.01;

		      // comment this line to implement adaptive time step
		      // kmax = 1;

		      dtshort = dtlong / (double) kmax;

		      /* store state of current point in U temporarily*/
		      for (m = 1; m <= num_states; m++)
            U[m] = u[n][m];

          /* integrate ODEs using Rush and Larsen scheme */
		      for (k = 1; k <= kmax; k++)
	    	    {
            if ((celltype[n] == 1) && (D[n] >= 0.025))
              dV = dtshort * calculate_TP06_current_OpSplit( U, dtshort, lookup, celltype[n], stimCurrent );
            else
              dV = 0.0;


 	    	    U[V] = U[V] - dV;
	    	    }

		      /* update state u with new values stored in U */
		      for (m = 1; m <= num_states; m++)
	   	 	    u[n][m] = U[m];

        }
/* end of step 2 */

/* step 3 */

      /* calculate diffusion */
      for (n = 1; n <= N; n++)
        {
        old_Vm[n] = new_Vm[n];
        dummy1 = u[n][V];
        dummy2 = dummy1 + half_dtlong * diffusion_2D( u, nneighb, n, N, D, dx2 );
        new_Vm[n] = dummy2;
        }

      /* update */
      for (n = 1; n <= N; n++)
        {
        dummy1 = new_Vm[n];
        u[n][V] = dummy1;
        }

      /* calculate diffusion */
      for (n = 1; n <= N; n++)
        {
        dummy1 = u[n][V];
        dummy2 = dummy1 + half_dtlong * diffusion_2D( u, nneighb, n, N, D, dx2 );
        new_Vm[n] = dummy2;
        dVdt[n] = dummy2 - old_Vm[n];
        }

/* end of step  3*/

/* detect upstrokes and downstrokes */

      if (time >= lastS1)
        {
        for (n = 1; n <= N; n++)
          {
          if ((new_Vm[n] > threshold) && (old_Vm[n] <= threshold) && (upStrokeTime[n][beat[n]] < 0) && (D[n] >= 0.025))
            {
            upStrokeTime[n][beat[n]] = time;
            }

    // only detect downstroke if upstroke has been detected first
          if ((new_Vm[n] < threshold) && (old_Vm[n] >= threshold) && (upStrokeTime[n][beat[n]] > 0) && (D[n] >= 0.025))
            {
            downStrokeTime[n][beat[n]] = time;

            if (beat[n] < numS1Beats)
              {
              beat[n]++;
              }
            }
          }
        }

/* output electrogram data every 1 ms */
	  if (modf(time/1.0, &timems) == 0.0)
		  {
      printf("time %f ms, writing electrograms to file\n",timems);

/* and write electrograms to eg file */
      fprintf(egPtr, "%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n", new_Vm[402],new_Vm[30100],new_Vm[30200],new_Vm[45150],new_Vm[60100],new_Vm[600]);
      printf("%4.2f %4.2f %4.2f %4.2f %4.2f\n",
        new_Vm[30100],new_Vm[30200],new_Vm[45150],new_Vm[60100],new_Vm[60200]);
      }

/* output stf file every 10 ms */
    if (modf(time/10.0, &timems) < 0.0001)
      {
      printf("time %f ms, writing stffile\n",time);
      stfout_2D( u, geom, stfcount*10, nrows, ncols );
      stfcount++;
      }

/* write to checkpoint file if needed */
      if (CHKPT_WRITE)
        {
    	  if (t == CHKPT_WRITE_TIME)
    	    {
    		dummy = checkpoint_write( u, time, t, CHKPT_WRITE_TIME, N );
    		exit(1);
    	    }
        }

  }
  printf("leaving main loop\n");

  /* save upstroke and downstroke data to files */
  for (n = 1; n <= N; n++)
    timing[n] = upStrokeTime[n][1];
  sprintf( outputFile,"%supStrokeTimeS1.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = downStrokeTime[n][1];
  sprintf( outputFile,"%sdownStrokeTimeS1.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = upStrokeTime[n][2];
  sprintf( outputFile,"%supStrokeTimeS2.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = downStrokeTime[n][2];
  sprintf( outputFile,"%sdownStrokeTimeS2.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = upStrokeTime[n][3];
  sprintf( outputFile,"%supStrokeTimeS3.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = downStrokeTime[n][3];
  sprintf( outputFile,"%sdownStrokeTimeS3.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = upStrokeTime[n][4];
  sprintf( outputFile,"%supStrokeTimeS4.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);

  for (n = 1; n <= N; n++)
    timing[n] = downStrokeTime[n][4];
  sprintf( outputFile,"%sdownStrokeTimeS4.stf",OUTPUTFILEROOT);
  writeData(outputFile, timing, geom, nrows, ncols);
  
  sprintf( outputFile,"%sdiffusion.stf",OUTPUTFILEROOT);
  writeData(outputFile, D, geom, nrows, ncols);

  fclose(egPtr);

  /* end of main loop */

  /* Free memory */
  free_fmatrix(lookup,0,num_lookup,0,voltage_steps);
  free_fmatrix(u,1,N,1,num_states);
  free_imatrix(geom, 1, nrows, 1, ncols);
  free_imatrix(nneighb, 1, RC, 1, 8);
  free_fvector(dVdt, 1, N );
  free_fvector(new_Vm, 1, N );
  free_fvector(old_Vm, 1, N );
  free_fvector(U, 1, RC);
  free_fvector(params, 1, num_params);
  free_ivector(celltype, 1, N);

  free_fmatrix(upStrokeTime, 1, N, 1, numS1Beats);
  free_fmatrix(downStrokeTime, 1, N, 1, numS1Beats);
  free_ivector(beat, 1, N);
  free_ivector(rowList, 1, N);
  free_ivector(colList, 1, N);

}

void writeData(char *fname, double *dataToWrite, int **geomarray, int nrows, int ncols)
{
    int lay, row, col;
    int index, outint;
    double outdouble;
    FILE *stf_file;

    stf_file = fopen( fname, "w" );

    fprintf(stf_file, "NAME Vm\n");
    fprintf(stf_file, "RANK 2\n");
    fprintf(stf_file, "DIMENSIONS %d %d\n",ncols, nrows);
    fprintf(stf_file, "BOUNDS %d %d %d %d\n", 0, ncols-1, 0, nrows-1);
    fprintf(stf_file, "SCALAR\n");
    fprintf(stf_file, "DATA\n");

    for (row = 1;row <= nrows; row+=1)
        {
        for (col = 1; col <= ncols; col+=1)
            {
            index = (geomarray[row][col] > 0)?geomarray[row][col]:0;
            if (index > 0)
                outdouble = dataToWrite[index];
            else
                outdouble = 0.0;
            fprintf(stf_file, "%6.4f ", outdouble);
            }
        fprintf( stf_file, "\n");
        }
    fclose(stf_file);
}
