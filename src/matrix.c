#include "matrix.h"
#include <stdio.h>

/*
 * Z = X^T Y
 *
 * X: m by n
 * Y: m by p
 * Z: n by p
 *
 * row major ordering
 *
 * This is a bit wasteful for Z = X^T X due to symmetry
 */
void crossprod(double *x, double *y, double *z, int m, int n, int p)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i * p + j] = x[k * n + i] * y[k * p + j];
	 for(k = 1 ; k < m ; k++)
	    z[i * p + j] += x[k * n + i] * y[k * p + j];
      }
   }
}

void cov(double *x, double *S, double *mean, int n, int p)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 S[i * p + j] = (S[k * n + i] - mean[j]) 
		  * (x[k * p + j] - mean[j]);
	 for(k = 1 ; k < m ; k++)
	    S[i * p + j] += (S[k * n + i] - mean[j]) 
		  * (S[k * p + j] - mean[j]);
      }
   }

   for(i = 0 ; i < p2 ; i++)
      S[i] *= d;
}

/* 
 * Invert 2x2 matrix x and put it in y, row-major matrix ordering:
 *  0 1
 *  2 3
 */
void invert2x2(double *y, double *x)
{
   double s = x[0] * x[3] - x[1] * x[2];

   y[0] =  x[3] / s;
   y[1] = -x[1] / s;
   y[2] = -x[2] / s;
   y[3] =  x[0] / s;
}

/*
 * Weighted cross product
 * Z = X^T W Y
 *
 * W: a diagonal m by m matrix
 * X: m by n
 * Y: m by p
 * Z: n by p
 *
 * Note that w is NOT actually a matrix, it's an array of length m
 */
void wcrossprod(double *x, double *y, double *w, double *z,
      int m, int n, int p)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i * p + j] = x[k * n + i] * y[k * p + j] * w[k];
	 for(k = 1 ; k < m ; k++)
	    z[i * p + j] += x[k * n + i] * y[k * p + j] * w[k];
      }
   }
}

/*
 * Square-Matrix-vector product
 *
 * z = X y
 *
 * z: m by 1
 * X: m by m
 * y: m by 1
 *
 */
void sqmvprod(double *x, double *y, double *z, int m)
{
   int i, k;
      
   for(i = 0 ; i < m ; i++)
   {
      k = 0;
      z[i] = x[k * m + i] * y[k];
      for(k = 1 ; k < m ; k++)
	 z[i] += x[k * m + i] * y[k];
   }
}

/* 
 * Copy the n by p matrix x to the n by m matrix y, ignoring p-m columns
 *
 */
void copyshrink(double *x, double *y, int n, int p, int *active, int m)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
   {
      k = 0;
      for(j = 0 ; j < p ; j++)
      {
	 /*printf("i:%d j:%d k:%d m:%d active[%d]:%d\n", i, j, k, m, j, active[j]);
	 fflush(stdout);*/
	 if(active[j])
	 {
	    y[i * m + k] = x[i * p + j];
	    k++;
	 }
      }
   }
}

/* 
 * Copy the n by p matrix x to the n by m matrix y
 *
 * from: start (inclusive)
 * to: end (exclusive)
 *
 */
void copyshrinkrange(double *x, double *y, int n, int p, int from, int to)
{
   int i, j, k, m = to - from;

   for(i = 0 ; i < n ; i++)
   {
      k = 0;
      for(j = 0 ; j < p ; j++)
      {
	 if(j >= from && j < to)
	 {
	    y[i * m + k] = x[i * p + j];
	    k++;
	 }
      }
   }
}

void printmatrix(double *x, int n, int p)
{
   int i, j;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p - 1; j++)
	 printf("%.3f ", x[i * p + j]);
      printf("%.3f\n", x[i * p + j]);
   }
}

