#include "pcor.h"

int writematrix(double **x, int n, int p, char* file)
{
   int i, j;
   FILE* out = fopen(file, "wt");
   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 fprintf(out, "%.20f", x[i][j]);
	 if(j < p - 1)
	    fprintf(out, ",");
      }
      fprintf(out, "\n");
   }

   fflush(out);
   fclose(out);
   return SUCCESS;
}

double sgd_matrix(double *beta, double *mean, double *sd,
      dloss dloss_func, loss_pt loss_pt_func, predict_pt predict_pt_func,
      double **x, int whichy, int n, int p)
{
   int i, j, k;
   int epoch = 0, maxepoch = 1000;
   double *grad = malloc((p - 1) * sizeof(double));
   double stepsize = 1e-5;
   double lambda1 = 1e-6, lambda2 = 0;
   double ptloss = 0;
   double *x2 = malloc((p - 1) * sizeof(double)); /* p-1 coefs  */
   double loss = 0, prevloss = 0;
   double y;
   
   while(epoch <= maxepoch)
   {
      loss = 0;
      for(i = 0 ; i < n ; i++)
      {
	 y = x[i][whichy];

	 k = 0;
	 for(j = 0 ; j < p ; j++)
	 {
	    if(j != whichy)
	    {
	       x2[k] = x[i][j];
	       k++;
	    }
	 }

	 ptloss = loss_pt_func(x2, beta, y, p - 1); 
	 /*yhat = predict_pt_func(&sm, beta, g->mean, g->sd, g->p + 1);*/
	 
	 dloss_func(x2, beta, y, p - 1, grad);
	 loss += ptloss;

	 /* Update weights */
	 for(j = 0 ; j < p - 1; j++)
	 {
	    beta[j] -= stepsize * (grad[j] 
	       + lambda1 * sign(beta[j]) 
	       + lambda2 * beta[j] * beta[j]);
	 }
      }
      printf("Epoch %d loss=%.6f\n", epoch, loss / n);
      if(epoch > 1 && fabs(prevloss - loss) <= 1e-9)
	 break;
      prevloss = loss;
      epoch++;
   }


   free(grad);
   free(x2);
   return 0;
}

void scale(double **x, double *mean, double *sd, int n, int p)
{
   int i, j;
   double delta;
 
   /* sd is really the sum of squares, not the SD, but we
    * use the same variable to save memory */

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
         if(i == 0)
	    mean[j] = sd[j] = 0;

         delta = x[i][j] - mean[j];
         mean[j] += delta / (i + 1);
         sd[j] += delta * (x[i][j] - mean[j]);
      }
   }

   for(j = 0 ; j < p ; j++)
      sd[j] = sqrt(sd[j] / (n - 1));
}


/*
 * Convert a p * (p - 1) matrix of regression coefficients (including the intercept)
 * to a p * p matrix of partial correlations with 1 on the diagonal
 */
void reg2pcor(double **beta, double **r, int p)
{
   int i, j, k;
   double **tmp = malloc(sizeof(double*) * p);

   for(i = 0 ; i < p ; i++)
   {
      tmp[i] = calloc(p, sizeof(double));
      k = 0;
      for(j = 0 ; j < p ; j++)
      {
	 if(j != i)
	 {
	    tmp[i][j] = beta[i][k];
	    k++;
	 }
      }
   }

   writematrix(tmp, p, p, "tmp.csv");

   /* convert regression coefs to partial correlation, shrink to zero if signs
    * don't agree */
   for(i = 0 ; i < p ; i++)
   {
      r[i][i] = 1;
      for(j = 0 ; j < i ; j++)
	 r[i][j] = r[j][i] = sign(tmp[i][j]) * 
	    sqrt(fmax(tmp[i][j] * tmp[j][i], 0));
   }

   for(i = 0 ; i < p ; i++)
      free(tmp[i]);
   free(tmp);
}

int main()
{
   int i, j;
   int n = 100;
   int p = 500;
   /* long seed = time(NULL); */
   long seed = 123;
   double **x = malloc(sizeof(double*) * n);
   double **x2 = malloc(sizeof(double*) * n);
   double **beta = malloc(sizeof(double*) * p);
   double *mean = malloc(sizeof(double) * p);
   double *sd = malloc(sizeof(double) * p);
   double **pcor = malloc(sizeof(double*) * p);
   gmatrix g;

   gmatrix_init(&g, TRUE, FALSE, NULL, x, NULL, n, p);
   
   srand48(seed);

   for(i = 0 ; i < n ; i++)
   {
      x[i] = malloc(sizeof(double) * p);
      for(j = 0 ; j < p ; j++)
	 x[i][j] = drand48();
   }

   /* scale so we don't have to worry about intercept */
   scale(x, mean, sd, n, p);

   for(i = 0 ; i < n ; i++)
   {
      x2[i] = malloc(sizeof(double) * p);
      for(j = 0 ; j < p ; j++)
	 x2[i][j] = (x[i][j] - mean[j]) / sd[j];
   }

   for(i = 0 ; i < p ; i++)
   {
      beta[i] = calloc(p-1, sizeof(double));
      pcor[i] = calloc(p, sizeof(double));
   }

   for(i = 0 ; i < p ; i++)
   {
      printf("%d\n", i);
      sgd_matrix(beta[i], mean, sd,
	    l2dloss, l2loss_pt, predict_l2loss_pt, x2, i, n, p);
   }

   writematrix(x2, n, p, "x.csv");
   writematrix(beta, p, p-1, "beta_pcor.csv");
   reg2pcor(beta, pcor, p);
   writematrix(pcor, p, p, "pcor.csv");

   return EXIT_SUCCESS;
}

