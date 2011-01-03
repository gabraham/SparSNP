#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "thin.h"
#include "matrix.h"

void cov(double *x, double *S, int n, int p)
{
   int i, m = n * p;
   double d = 1.0 / (n - 1.0);

   crossprod(x, x, S, n, p, p);
   for(i = 0 ; i < m ; i++)
      S[i] *= d;
}

/* Converts a p by p covariance matrix S to a p by p correlation matrix P
 */
void cov2cor(double *S, double *P, int p)
{
   double z;
   int i, j;

   for(i = 1 ; i < p ; i++)
   {
      for(j = 0 ; j < i ; j++)
      {
	 z = S[i * p + j] / sqrt(S[i * p + i] * S[j * p + j]);
	 P[i * p + j] = P[j * p + i] = z;
      }
   }

   for(i = 0 ; i < p ; i++)
      P[i * p + i] = 1.0;
}

/* Remove SNPs based on correlation in a sliding window of size THIN_WINDOW_SIZE
 */
int thin(double *x, int n, int p, int *active,
      double cormax, int windowsize, int stepsize)
{
   int j, i, k;
   double *S = NULL,
	  *P = NULL,
	  *xwin = NULL;
   int ws = THIN_WINDOW_SIZE,
       ws2 = ws * ws;


   j = 0;
   while(j < p && ws > 1)
   {
      MALLOCTEST(xwin, sizeof(double) * ws * n);
      copyshrink(x, xwin, n, p, active, ws);

      MALLOCTEST(S, sizeof(double) * ws2);
      MALLOCTEST(P, sizeof(double) * ws2);

      /* correlation in the window */
      cov(x, S, n, p);
      FREENULL(xwin);
      cov2cor(S, P, ws);
      FREENULL(S);

      /* filter pairs looking at lower triangular matrix */
      for(i = 1 ; i < ws ; i++)
      {
	 for(k = 0 ; k < i; k++)
	 {
	    /* remove the first member of each pair, if both are active */
	    if(fabs(P[i * ws + k]) > cormax && active[i + j] && active[k + j])
	       active[i + j] = FALSE;
	 }
      }

      FREENULL(P);

      /* beware of edge effect */
      ws = (windowsize < p - j) ? windowsize : p - j;
      j += THIN_STEP_SIZE;
   }

   return SUCCESS;
}

