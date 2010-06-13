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

/* variant of SGD for partial correlation */
double sgd_matrix(double *beta, double *mean, double *sd,
      dloss_pt dloss_pt_func, loss_pt loss_pt_func, predict_pt predict_pt_func,
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

	 /*ptloss = loss_pt_func(x2, beta, y, p - 1);  */
	 /*yhat = predict_pt_func(&sm, beta, g->mean, g->sd, g->p + 1);*/
	 
	 /*dloss_func(x2, beta, y, p - 1, grad);*/
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

/*void scale(double **x, double *mean, double *sd, int n, int p)
{
   int i, j;
   double delta;
 
   / sd is really the sum of squares, not the SD, but we
     use the same variable to save memory /

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
}*/


/*
 * Convert a p * p matrix of regression coefficients (including the intercept)
 * to a p * p matrix of partial correlations with 1 on the diagonal
 */
void reg2pcor(double **beta, double **r, int p)
{
   int i, j, k;
   double **tmp;
   
   /*tmp = malloc(sizeof(double*) * p);

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
   }*/

   /*writematrix(tmp, p, p, "tmp.csv");*/

   tmp = beta;

   /* convert regression coefs to partial correlation, shrink to zero if signs
    * don't agree */
   for(i = 0 ; i < p ; i++)
   {
      r[i][i] = 1;
      for(j = 0 ; j < i ; j++)
	 r[i][j] = r[j][i] = sign(tmp[i][j]) * 
	    sqrt(fmax(tmp[i][j] * tmp[j][i], 0));
   }

   /*for(i = 0 ; i < p ; i++)
      free(tmp[i]);
   free(tmp);*/
}

int main(int argc, char *argv[])
{
   int i, j;
   char *filename = NULL;
   int n = 0;
   int p = 0;
   long seed = 123;
   double stepsize = 1e-4;
   double lambda1 = 0, lambda2 = 0;
   double threshold = 1e-6;
   short verbose = TRUE;
   double trunc = 1e-9;
   double **beta;
   double **pcor;
   int *trainf = NULL;
   optim_gmatrix optim_gmatrix_func = gd_gmatrix;
   loss_pt loss_pt_func = NULL;
   predict_pt predict_pt_func = NULL;
   dloss_pt dloss_pt_func = NULL;
   predict_gmatrix predict_gmatrix_func = NULL;
   gmatrix g;
   dtype *tmp;
   int maxepochs = 100;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-optim"))
      {
	 i++;
	 if(strcmp2(argv[i], "sgd"))
	    optim_gmatrix_func = sgd_gmatrix;
	 else if(strcmp2(argv[i], "scd"))
	    optim_gmatrix_func = scd_gmatrix;
	 else if(strcmp2(argv[i], "gd"))
	    optim_gmatrix_func = gd_gmatrix;
      }
      else if(strcmp2(argv[i], "-f"))
      {
	 i++;
	 filename = argv[i];
      }
      else if(strcmp2(argv[i], "-n"))
      {
	 i++;
	 n = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-p"))
      {
	 i++;
	 p = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-epochs"))
      {
	 i++;
	 maxepochs = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-step"))
      {
	 i++;
	 stepsize = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1"))
      {
	 i++;
	 lambda1 = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l2"))
      {
	 i++;
	 lambda2 = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-thresh"))
      {
	 i++;
	 threshold = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-v"))
      {
	 verbose = TRUE;
      }
      else if(strcmp2(argv[i], "-vv"))
      {
	 verbose = 2;
      }
   }

   if(filename == NULL || n == 0 || p == 0)
   {
      printf("usage: pcor -f <filename> -n <#samples> -p <#variables>\n");
      return EXIT_FAILURE;
   }


   if(!gmatrix_init(&g, TRUE, TRUE, filename, NULL, NULL, n, p - 1))
      return EXIT_FAILURE;
   
   gmatrix_scale(&g);
   writevectorf("mean.csv", g.mean, p);
   writevectorf("sd.csv", g.sd, p);


   srand48(seed);

   MALLOCTEST2(trainf, sizeof(int) * g.n)
   for(i = 0 ; i < g.n ; i++)
      trainf[i] = TRUE;

   loss_pt_func = &l2loss_pt;
   dloss_pt_func = &l2dloss_pt;
   predict_gmatrix_func = &predict_l2loss_gmatrix;
   predict_pt_func = &predict_l2loss_pt;

   MALLOCTEST2(beta, p * sizeof(double*))
   MALLOCTEST2(pcor, p * sizeof(double*))
   
   for(j = 0 ; j < p ; j++)
   {
      CALLOCTEST2(beta[j], p, sizeof(double))
      CALLOCTEST2(pcor[j], p, sizeof(double))
   }

   MALLOCTEST2(tmp, n * sizeof(dtype))
   
   for(j = 0 ; j < p ; j++)
   {
      printf("%d\n", j);


      /* set current variable to 1 (intercept) */
      for(i = 0 ; i < n ; i++)
      {
	 /* remember this variable and use it as response in current round */
	 tmp[i] = g.x[i][j];

	 /* the intercept is the jth variable, ends up being the diagonal of
	  * the beta matrix
	  */
	 g.x[i][j] = 1.0;
      }
      g.y = tmp;

      optim_gmatrix_func(&g, dloss_pt_func, loss_pt_func, predict_pt_func,
	 stepsize, maxepochs, beta[j], lambda1, lambda2, threshold,
	 verbose, trainf, trunc);

      /* put variable back */
      for(i = 0 ; i < n ; i++)
	 g.x[i][j] = tmp[i];
   }


   writematrix(beta, p, p, "beta_pcor.csv");
   reg2pcor(beta, pcor, p);
   writematrix(pcor, p, p, "pcor.csv");

   for(i = 0 ; i < p ; i++)
   {
      free(beta[i]);
      free(pcor[i]);
   }

   free(beta);
   free(pcor);


   return EXIT_SUCCESS;
}

