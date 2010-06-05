#include "sgd.h"

/* Stochastic gradient descent */
double sgd_gmatrix(gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   loss_pt loss_pt_func,    /* loss for one sample */
   predict_pt predict_pt_func, /* prediction for one sample */
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc)
{
   int epoch = 1, i, j;
   double *grad;
   double prevloss = 0, loss = 0, bestloss = 1e9;
   double stepsize = maxstepsize;
   sample sm;
   double yhat;
   double trainacc = 0;
   double testacc = 0;
   double testloss = 0;
   double diff = 1e9;
   int ntests = 0;
   double ptloss = 0;
   double s, d, dp = 0;

   MALLOCTEST(grad, sizeof(double) * (g->p + 1))

   if(!sample_init(&sm, g->p))
      return FAILURE;

   while(epoch <= maxepoch)
   {
      loss = 0;
      testloss = 0;
      testacc = 0;
      trainacc = 0;
      ntests = 0;

      for(i = 0 ; i < g->n ; i++)
      { 
	 g->nextrow(g, &sm);

	 dp = dotprod(sm.x, beta, g->p + 1);
	 ptloss = loss_pt_func(dp, sm.y); 
	 yhat = predict_pt_func(dp);

	 /* train */
	 if(trainf[i])
	 {
	    dloss_pt_func(sm.x, dp, sm.y, g->p + 1, grad);
	    loss += ptloss;
	    trainacc += (double)((yhat >= 0.5) == (int)sm.y);

	    /* Update weights */
	    for(j = 0 ; j < g->p + 1; j++)
	    {
	       s = stepsize;
	       d = grad[j] + lambda1 * sign(beta[j]) 
		     + lambda2 * beta[j] * beta[j];
	       if(sign(beta[j]) != sign(d))
		  s = stepsize / 2.0;
	       beta[j] -= s * d;
	    }
	 }
	 /* test */
	 else
	 {
	    testloss += ptloss;
	    testacc += (double)((yhat >= 0.5) == (int)sm.y);
	    ntests++;
	 }
      }

      trainacc = trainacc / (g->n - ntests);

      if(ntests > 0)
      {
	 testacc = testacc / ntests;
	 testloss = testloss / ntests;
      }
      /*printf("total loss: %.10f over %d samples\n", loss, g->n - ntests);*/
      loss = loss / (g->n - ntests);
      if(bestloss > loss)
	 bestloss = loss;

      /* truncate small weights when lasso is active */
      if(lambda1 > 0)
         for(j = 0 ; j < g->p + 1; j++)
	    if(fabs(beta[j]) <= trunc)
	       beta[j] = 0;

      diff = prevloss - loss;

      if(verbose)
      {
	 printf("Epoch %d  training loss: %.5f diff: %.5f stepsize: %.15f\
 training accuracy: %.8f", epoch, loss, diff, stepsize, trainacc);
	 if(ntests > 0)
	    printf(" test accuracy: %.8f test loss: %.8f", testacc, testloss);
	 printf("\n");
      }
 
      if(epoch > 1 && diff < -threshold)
	 stepsize = fmax(stepsize / 2, 1e-20);
      else if(fabs(diff) <= threshold)
      {
	 if(verbose)
	    printf("Termination condition met\n");
	 break;
      }

      prevloss = loss;
      epoch++;
   }

   sample_free(&sm);
   free(grad);
   return loss;
}

void writeout(char* file, double **x, int *y, double n, double p)
{
   int i, j;
   FILE* out = fopen(file, "w+");

   for(i = 0 ; i < n ; i++)
   {
      fprintf(out, "%d,", y[i]);
      for(j = 0 ; j < p ; j++)
      {
	 if(j < p - 1)
	    fprintf(out, "%.20f,", x[i][j]);
	 else
	    fprintf(out, "%.20f\n", x[i][j]);
      }
   }

   fflush(out);
   fclose(out);
}

