#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "thin.h"
#include "matrix.h"

void cov(double *x, double *S, int n, int p)
{
   int i, p2 = p * p;
   double d = 1.0 / (n - 1.0);

   crossprod(x, x, S, n, p, p);
   for(i = 0 ; i < p2 ; i++)
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
 *
 * x is the original n by p data
 */
int thin(double *x, int n, int p, int *active,
      double cormax, int windowsize, int stepsize)
{
   int j, i, k;
   double *S = NULL,
	  *P = NULL,
	  *xwin = NULL;
   int ws = windowsize,
       ws2 = ws * ws;

   j = 0;
   while(j < p && ws > 1)
   {
      printf("%d\n", j);
      fflush(stdout);

      /* Consider SNPs in the window only */
      MALLOCTEST(xwin, sizeof(double) * ws * n);
      copyshrinkrange(x, xwin, n, p, j, j + ws);

      MALLOCTEST(S, sizeof(double) * ws2);
      MALLOCTEST(P, sizeof(double) * ws2);

      printf("ws: %d\tws2: %d\n", ws, ws2);
      fflush(stdout);
      /* Estimate correlation in the window */
      cov(xwin, S, n, ws);
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
      ws2 = ws * ws;
      j += stepsize;
   }

   return SUCCESS;
}

