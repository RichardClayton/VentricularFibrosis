
/***************************************************************

 nrutils.c

 Version       2.1

 Date          23-oct-2023

 Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)
 
 This file is part of VentricularFibrosis.

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

/* Numerical recipes routines */

/***************************************************************

 fvector
 NR routine to allocate space for a double vector with 
 subscript range v[nl..nh]

***************************************************************/

double *fvector( long nl, long nh )
{
   double *v;
   v = (double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in fvector()");
   return v-nl+NR_END;
}

/***************************************************************

Routine ivector - allocates space for an int vector with 
subscript range v[nl..nh]

***************************************************************/

int *ivector( long nl, long nh )
{
   int *v;

   v = (int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
   return v-nl+NR_END;
}


/**************************************************************

 imatrix
 NR routine to allocate space for a int matrix with 
 subscript range v[nrl..nrh][ncl..nch]

***************************************************************/

int **imatrix( long nrl, long nrh, long ncl, long nch )
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **mm;

  /* allocate pointers to rows */
  mm = (int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!mm) 
    nrerror("allocation failure 1 in matrix()");
  mm += NR_END;
  mm -= nrl;

  /* allocate rows and set pointers to them */
  mm[nrl] = (int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!mm[nrl])
    nrerror("allocation failure 2 in matrix()");
  mm[nrl] += NR_END;
  mm[nrl] -= ncl;

  for ( i = nrl+1; i<= nrh; i++ )
    mm[ i ] = mm[ i - 1 ] + ncol;

  /* return pointer to array of pointers to rows */

  return mm;
}


/**************************************************************

 fmatrix
 NR routine to allocate space for a int matrix with 
 subscript range v[nrl..nrh][ncl..nch]

***************************************************************/

double **fmatrix( long nrl, long nrh, long ncl, long nch )
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **mm;

  /* allocate pointers to rows */
  mm = (double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!mm) 
    nrerror("allocation failure 1 in matrix()");
  mm += NR_END;
  mm -= nrl;

  /* allocate rows and set pointers to them */
  mm[nrl] = (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!mm[nrl])
    nrerror("allocation failure 2 in matrix()");
  mm[nrl] += NR_END;
  mm[nrl] -= ncl;

  for ( i = nrl+1; i<= nrh; i++ )
    mm[ i ] = mm[ i - 1 ] + ncol;

  /* return pointer to array of pointers to rows */

  return mm;
}

/**************************************************************

 i3dmatrix
 NR routine to allocate space for an int 3d matrix with 
 subscript range v[nrl..nrh][ncl..nch][ndl..ndh]

***************************************************************/

int ***i3dmatrix( long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
  long i, jj, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
  int ***t;
  
/* allocate pointers to pointers to rows */
  t = (int ***) malloc((size_t)((nrow + NR_END) * sizeof(int **)));
  if (!t) nrerror("allocation failure 1 in i3dmatrix\n");
  t += NR_END;
  t -= nrl;

/* allocate pointers to rows and set pointers to them */
  t[nrl] = (int **) malloc((size_t)((nrow * ncol + NR_END)*sizeof(int *)));
  if (!t[nrl]) nrerror("allocation failure 2 in i3dmatrix\n");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

/* allocate rows and set pointers to them */
  t[nrl][ncl] = (int *) malloc((size_t)((nrow * ncol * ndep + NR_END)*sizeof(int)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in i3dmatrix\n");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(jj = ncl + 1;jj <= nch; jj++) t[nrl][jj] = t[nrl][jj-1] + ndep;
  for(i = nrl + 1;i <= nrh; i++)
    {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol * ndep;
    for(jj = ncl + 1;jj <= nch; jj++) t[i][jj] = t[i][jj-1] + ndep;
    }

/* return pointer to array of pointers to rows */

  return t;
}

/**************************************************************

 free_ivector
 NR routine to free double vector with subscript range
 v[nl..nh]

***************************************************************/

void free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/**************************************************************

 free_imatrix
 NR routine to free int matrix with subscript range
 v[nrl..nrh][ncl..nch]

***************************************************************/

void free_imatrix( int **mm, long nrl, long nrh, long ncl, long nch )
{
  free((FREE_ARG) (mm[nrl]+ncl-NR_END));
  free((FREE_ARG) (mm+nrl-NR_END));
}

/**************************************************************

 free_fvector
 NR routine to free double vector with subscript range
 v[nl..nh]

***************************************************************/

void free_fvector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/**************************************************************

 free_fmatrix
 NR routine to free double matrix with subscript range
 v[nrl..nrh][ncl..nch]

***************************************************************/

void free_fmatrix( double **mm, long nrl, long nrh, long ncl, long nch )
{
  free((FREE_ARG) (mm[nrl]+ncl-NR_END));
  free((FREE_ARG) (mm+nrl-NR_END));
}

/**************************************************************

 free_i3dmatrix
 NR routine to free int 3dmatrix with subscript range
 v[nrl..nrh][ncl..nch][ndl..ndh]

***************************************************************/

void free_i3dmatrix( int ***mm, long nrl, long nrh, long ncl, long nch, 
                                                          long ndl, long ndh )
{
  free((FREE_ARG) (mm[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (mm[nrl]+ncl-NR_END));
  free((FREE_ARG) (mm+nrl-NR_END));
}

/**************************************************************

 Error handler

***************************************************************/

void nrerror(char error_text[])
{
  fprintf(stderr,"error: ");
  fprintf(stderr,"%s\n", error_text);
  exit(1);
}

