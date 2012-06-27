#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "thin.h"
#include "matrix.h"

/* Remove SNPs based on correlation in a sliding window.
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

   j = 0;
   while(j < p && ws > 1)
   {
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
	    /* remove the first member of each pair, if both are active */
	    if(fabs(P[i * ws + k]) > cormax && active[i + j] && active[k + j])
	    {
	       active[i + j] = FALSE;
	       /*printf("var %d filtered, P[%d,%d]: %.5f\n", i + j, i, k, P[i * ws + k]);*/
	    }
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

