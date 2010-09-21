#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"
#include "scale.h"

/*
 * Converts a matrix encoded as chars (one byte per entry), scales it, and
 * saves it as a matrix of doubles (8 bytes per entry).
 *
 * Expects data in column major ordering
 */
int scale(gmatrix *g)
{
   int i, j, p1 = g->p + 1;
   sample sm;
   double *tmp = NULL, delta;

   if(!sample_init(&sm, g->ncurr))
      return FAILURE;

   if(!g->mean)
      MALLOCTEST(g->mean, sizeof(double) * p1);
   if(!g->sd)
      MALLOCTEST(g->sd, sizeof(double) * p1);
   g->mean[0] = 0;
   g->sd[0] = 1;

   MALLOCTEST(tmp, sizeof(double) * g->ncurr);
   
   /* read intercept and ignore it*/
   g->nextcol(g, &sm, 0);

   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d of %d", j, p1);
      g->nextcol(g, &sm, j);
      g->mean[j] = g->sd[j] = 0;
      for(i = 0 ; i < g->ncurr ; i++)
      {
	 delta = sm.x[i] - g->mean[j];
	 g->mean[j] += delta / (i + 1);
	 g->sd[j] += delta * (sm.x[i] - g->mean[j]);
      }

      g->sd[j] = sqrt(g->sd[j] / (g->ncurr - 1));
      printf("\r");
   }
   printf("\n");

   free(tmp);
   sample_free(&sm);

   return SUCCESS;
}

/* unscales the vector of coefficients */
void unscale_beta(double *beta, double *mean, double *sd, int p)
{
   int j;

   beta[0] = 1.0;

   for(j = 1 ; j < p ; j++)
   {
      if(sd[j] == 0)
	 beta[j] = beta[j] + mean[j];
      else
	 beta[j] = beta[j] * sd[j] + mean[j];
   }
}

int writescale(char* filename, double *mean, double *sd, int p)
{
   FILE *out;
   FOPENTEST(out, filename, "wb")
   FWRITETEST(mean, sizeof(double), p, out);
   FWRITETEST(sd, sizeof(double), p, out);
   fclose(out);
   return SUCCESS;
}

int readscale(char* filename, double *mean, double *sd, int p)
{
   FILE *in;
   FOPENTEST(in, filename, "rb")
   FREADTEST(mean, sizeof(double), p, in);
   FREADTEST(sd, sizeof(double), p, in);
   fclose(in);
   return SUCCESS;
}

int main(int argc, char* argv[])
{
   int i, n = 0, p = 0, len, nfolds = 1, k;
   char *filename_bin = NULL,
        *filename_scale = "scale.bin",
	*filename_beta = NULL,
	*filename_beta_out = NULL,
	*filename_folds_ind = NULL;
   short doscale = TRUE,
	 encoded = TRUE,
	 binformat = BINFORMAT_BIN;
   gmatrix g;
   char tmp[100];

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-unscale"))
	 doscale = FALSE;
      if(strcmp2(argv[i], "-bin"))
      {
	 i++;
	 filename_bin = argv[i];
      }
      else if(strcmp2(argv[i], "-scale"))
      {
	 i++;
	 filename_scale = argv[i];
      }
      else if(strcmp2(argv[i], "-betafile"))
      {
	 i++;
	 filename_beta = argv[i];
	 len = strlen(filename_beta) + 1 + 5;
	 MALLOCTEST(filename_beta_out, len);
	 snprintf(filename_beta_out, len, "%s.unsc", filename_beta);
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
      else if(strcmp2(argv[i], "-notencoded"))
	 encoded = FALSE;
      else if(strcmp2(argv[i], "-plink"))
	 binformat = BINFORMAT_PLINK;
      else if(strcmp2(argv[i], "-foldind"))
      {
	 i++;
	 filename_folds_ind = argv[i];
      }
   }

   if((doscale && (filename_bin == NULL || n == 0 || p == 0))
      || (!doscale && 
	    (filename_scale == NULL || p == 0 || filename_beta == NULL)))
   {
      printf("scale: [-unscale] -bin <filein> [-scale <fileout>] \
[-betafile <betafile>] [-notencoded] [-plink] -n #n -p #p \
[-foldind <folds ind file>]\n");
      return EXIT_FAILURE;
   }

   if(doscale)
   {
      if(!gmatrix_init(&g, filename_bin, n, p,
	    NULL, YFORMAT01, MODEL_LINEAR,
	    encoded, binformat, filename_folds_ind, nfolds, MODE_TRAIN))
	 return EXIT_FAILURE;

      if(filename_folds_ind)
      {
	 for(k = 0 ; k < g.nfolds ; k++)
	 {
	    gmatrix_set_fold(&g, k);

	    if(!scale(&g))
	       return EXIT_FAILURE;

	    len = strlen(filename_scale) + 1 + 3;
	    snprintf(tmp, len, "%s.%02d", filename_scale, k);
	    if(!writescale(tmp, g.mean, g.sd, p + 1))
	       return EXIT_FAILURE;
	 }
      }
      else
      {
	 if(!scale(&g))
	    return EXIT_FAILURE;
	 if(!writescale(filename_scale, g.mean, g.sd, p + 1))
	    return EXIT_FAILURE;
      }

      gmatrix_free(&g);
   }
   else /* unscale beta coefficients, i.e., put on original scale of
	   unstandardised data */
   {
      g.p = p;
      MALLOCTEST(g.mean, sizeof(double) * (p + 1));
      MALLOCTEST(g.sd, sizeof(double) * (p + 1));
      MALLOCTEST(g.beta, sizeof(double) * (p + 1));

      if(!readscale(filename_scale, g.mean, g.sd, p + 1))
	 return EXIT_FAILURE;

      if(!load_beta(g.beta, filename_beta, g.p + 1))
	 return FAILURE;
 
      unscale_beta(g.beta, g.mean, g.sd, g.p + 1);

      if(!writevectorf(filename_beta_out, g.beta, g.p + 1))
	 return FAILURE;

      free(g.mean);
      free(g.sd);
      free(g.beta);
   }
   
   if(filename_beta_out)
      free(filename_beta_out);

   return EXIT_SUCCESS;
}

