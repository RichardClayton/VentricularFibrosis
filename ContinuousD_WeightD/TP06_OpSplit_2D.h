/***************************************************************

TP06_OpSplit_2D.h

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

/* Header file for TP06 in 2D with isotropic diffusion */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define ROWS            400      //340
#define COLUMNS         400
#define NUM_STATES      20      /* number of states in u array */
#define NUM_ITERATIONS  30000   /* maximum duration of simulation (60000 is 6s) */
//#define NUM_POINTS      100   /* Length of cable */
#define DT              0.10    /* maximum timestep (ms) */
#define DX              0.25    /* 0.25 mm : Assume dx = dy */
#define DIFFUSION       0.2     /* mm2/ms :isotropic diffusion */
#define VMOFFSET        1001.0  /* offset for lookup table */
#define GAIN            10.0    /* gain for lookup table */
#define NUM_LOOKUP      40      /* number of variables in lookup table*/
#define NUM_PARAMS      50      /* number of cell model parameters */
#define VOLTAGE_STEPS   2001    /* number of voltage steps in lookup table */
#define CM              1.0     /* Membrane capacitance (uF/cm2) */

#define S1BCL           400.0   /* basic cycle length for pacing (ms) */
#define NUMS1BEATS      3       /* number of S1 beats */
#define NUMS2BEATS      4       /* number of extra beats delivered at ERP */
//define S1S2            360.0   /* S1S2 interval (ms) 500 ms for SR, 400 for cAF*/
//#define DECREMENT       20.0    /* decrement for S1 beats */

/* Macros for numerical recipes routines */
#define FREE_ARG       	char*
#define NR_END         	1

/* diffusion file */
#define DIFFUSIONFILE   "DiffusionCoefficient.txt"

/* Macros for output */
#define OUTPUTFILEROOT  "TP06_2D_"
#define STFFILEROOT     "STFfiles/TP06_2D_"

/* checkpointing */
/* read/write times should be multiples  f 5 ms */
#define CHKPTROOT 		    "checkpoint"
#define CHKPT_READ	        0
#define CHKPT_READ_TIME     200000 // 4000 ms
#define CHKPT_WRITE	        0
#define CHKPT_WRITE_TIME    200000 // 4000 ms

/* forward declaration of all functions used */

/* PDE solver */
int initialise_geometry_2D( int **geom, int nrows, int ncols, int **nneighb, double *D );
void initialise_variables_2D( double **u, int N );
void initialise_spiral_2D( double **u, int N, int ny, int nx );
//void initialise_diffusion_2D( double *D, int nrows, int ncols );

int create_TP06_lookup_OpSplit_2D( double **lookup );
double calculate_TP06_current_OpSplit( double *U, double dt, double **lookup, int celltype, double stimCurrent );
double diffusion_2D_modD( double **u, int **nneighb, int n, int N, double *D, double dx2 );
int free_arrays( double **u, int N, int num_parameters, double **lookup, int voltage_steps, int num_lookup );

/* checkpointing */
int checkpoint_write( double **u, double time, int t, int count, int N );
int checkpoint_read( double **u, double *time, int *t, int count, int N);

/* Numerical recipes routines */
double *fvector( long nl, long nh );
int *ivector( long nl, long nh );
int **imatrix( long nrl, long nrh, long ncl, long nch );
double **fmatrix( long nrl, long nrh, long ncl, long nch );
int ***i3dmatrix( long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector( int *m, long nl, long nh );
void free_imatrix( int **m, long nrl, long nrh, long ncl, long nch );
void free_fvector( double *m, long nl, long nh );
void free_fmatrix( double **m, long nrl, long nrh, long ncl, long nch );
void free_i3dmatrix( int ***m, long nrl, long nrh, long ncl, long nch, long ndl, long ndh );
void nrerror( char error_text[] );

int stfout_2D( double **u, int **geomarray, int stfcount, int nx, int ny );

/* RGB file output */
//int write_rgb( int nx, int ny, int N, double **u, int **geom, char *fname, int I, double iMax, double iMin );
//void putbyte( FILE *outfile, unsigned char value );
//void putshort( FILE *outfile, unsigned short value );
//static int putlong( FILE *outfile, unsigned long value );

/*** END ***/
