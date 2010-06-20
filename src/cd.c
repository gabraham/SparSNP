#include "sgd.h"

/* coordinate descent */
double cd_gmatrix(gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   d2loss_pt d2loss_pt_func,        /* 2nd deriv */
   d2loss_pt_j d2loss_pt_j_func,        /* 2nd deriv wrt beta_j */
   loss_pt loss_pt_func,    /* loss for one sample */
   predict_pt predict_pt_func, /* prediction for one sample */
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc)
{
   int i, j, k;
   int epoch = 1;
   double d;
   double loss = 0;
   /*double *lp;*/
   double y2 = 0;
   double d2 = 0;
   double tmp;
   short *converged = NULL;
   int numconverged = 0;
   double relerr;

   if(!g->inmemory)
   {
      fprintf(stderr, "cd_gmatrix doesn't support disk based gmatrix yet");
      return FAILURE;
   }

   /*sample_init(&sm, g->n);*/

   CALLOCTEST(converged, g->p + 1, sizeof(short));

   while(epoch <= maxepoch && numconverged < g->p + 1)
   {
      printf("epoch %d\n", epoch);
      loss = 0;
      /*CALLOCTEST(lp, g->n, sizeof(double));*/

      for(j = 0 ; j < g->p + 1; j++)
      {
	 printf("%d", j);
	 fflush(stdout);
	 if(converged[j])
	    continue;

	 d = 0;
	 d2 = 0;

	 for(i = 0 ; i < g->n ; i++)
	 {
	    /*y2 = 0;
	    for(k = 0 ; k < g->p + 1 ; k++)
	       y2 += g->x[i][k] * beta[k];

	    d += g->x[i][j] * (g->y[i] - y2);*/
	    
	    /* For linear regression, the 2nd
	     * derivative wrt beta_j is \sum_{i=1}^n x_{ij}^2,
	     * which is always 1 / (n - 1) for standardised inputs
	     */
	    /*d2 += pow(g->x[i][j], 2);*/

	 }

	 /* TODO: don't penalise intercept */
	 tmp = soft_threshold(beta[j] + d / d2, lambda1) / (1 + lambda2);
	 /*tmp = beta[j] + d / d2;*/

	 if(epoch > 1)
	 {
	    relerr = fabs(beta[j] - tmp) / (fabs(beta[j]) + fabs(tmp));
	    if(relerr < threshold)
	    {
	       converged[j] = TRUE;
	       numconverged++;
	    }
	 }
	 beta[j] = tmp;

	 printf("\r");
      }

      /*for(i = 0 ; i < g->n ; i++)
      {
	 for(j = 0 ; j < g->p + 1 ; j++)
	    lp[i] += g->x[i][j] * beta[j];

	 loss += loss_pt_func(lp[i], g->y[i]) / g->n;
      }

      printf("Epoch %d  training loss: %.5f\n", epoch, loss);

      free(lp); */
      printf("Epoch %d done\n", epoch);
      epoch++;
   }

   free(converged);

   return SUCCESS;
}

