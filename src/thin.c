#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "thin.h"
#include "matrix.h"

/*int cov(double *x, double *S, int n, int p)
{
   int i, j, p2 = p * p;
   double d = 1.0 / (n - 1.0);
   double *mean = NULL;

   CALLOCTEST(mean, p, sizeof(double));
   for(j = 0 ; j < p ; j++)
   {
      for(i = 0 ; i < n ; i++)
	 mean[j] +=  x[i * p + j];
      mean[j] /= n;
   }

   crossprod(x, x, S, n, p, p);
   for(i = 0 ; i < p2 ; i++)
      S[i] *= d;

   FREENULL(mean);
}*/

/* Remove SNPs based on correlation in a sliding window of size THIN_WINDOW_SIZE
 *
 * x is the original n by p data
 *
 * The selected SNPs are indicated by the array "active"
 */
int thin(double *x, int n, int p, int *active,
      double cormax, int windowsize, int stepsize)
{
   int j, i, k;
   double *S = NULL,
	  *P = NULL,
	  *xwin = NULL;
   int ws = (windowsize < p) ? windowsize : p,
       ws2 = ws * ws;

   stepsize = (stepsize > windowsize) ? windowsize : stepsize;

   printf(">>>>p=%d\n", p); fflush(stdout);
   j = 0;
   while(j < p && ws > 1)
   {
      printf("ws=%d\tstepsize=%d\n", ws, stepsize); fflush(stdout);
      /* Consider SNPs in the window only */
      MALLOCTEST(xwin, sizeof(double) * ws * n);
      copyshrinkrange(x, xwin, n, p, j, j + ws);

      MALLOCTEST(S, sizeof(double) * ws2);
      MALLOCTEST(P, sizeof(double) * ws2);

      /* Estimate correlation in the window */
      cov(xwin, S, n, ws);
      cov2cor(S, P, ws);

      FREENULL(xwin);
      FREENULL(S);

      /* filter pairs looking at lower triangular matrix */
      for(i = 1 ; i < ws ; i++)
      {
	 for(k = 0 ; k < i; k++)
	 {
	    printf("i:%d j:%d k:%d ws:%d\n", i, j, k, ws); fflush(stdout);
	    /* remove the first member of each pair, if both are active */
	    if(fabs(P[i * ws + k]) > cormax && active[i + j] && active[k + j])
	       active[i + j] = FALSE;
	 }
      }

      FREENULL(P);

      /* beware of edge effect */
      stepsize = (stepsize < p - j) ? stepsize : p - j;
      j += stepsize;
      ws = (windowsize < p - j) ? windowsize : p - j;
      ws2 = ws * ws;
   }

   return SUCCESS;
}

