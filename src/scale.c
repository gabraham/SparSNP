/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include "common.h"
#include "cd.h"
#include "util.h"

/*
 * Converts a matrix encoded as chars (one byte per entry),
 * scales it, and saves it as a matrix of doubles (8 bytes
 * per entry).
 *
 * Expects data in column major ordering
 */
int scale(gmatrix *g)
{
   int i, j, p1 = g->p + 1, n = g->ncurr, ngood = 0;
   sample sm;
   double *tmp = NULL, delta;

   if(!sample_init(&sm))
      return FAILURE;

   if(!g->mean)
      MALLOCTEST(g->mean, sizeof(double) * p1);
   if(!g->sd)
      MALLOCTEST(g->sd, sizeof(double) * p1);
   g->mean[0] = 0;
   g->sd[0] = 1;

   MALLOCTEST(tmp, sizeof(double) * n);
   
   /* read intercept and ignore it*/
   g->nextcol(g, &sm, 0, NA_ACTION_ZERO);

   for(j = 1 ; j < p1 ; j++)
   {
      /*printf("%d of %d", j, p1);*/
      g->nextcol(g, &sm, j, NA_ACTION_ZERO);
      ngood = g->mean[j] = g->sd[j] = 0;
      for(i = 0 ; i < n ; i++)
      {
	 /* skip missing observations */
	 if(sm.x[i] != X_LEVEL_NA)
	 {
	    delta = sm.x[i] - g->mean[j];
	    g->mean[j] += delta / (i + 1);
	    g->sd[j] += delta * (sm.x[i] - g->mean[j]);
	    ngood++;
	 }
      }

      g->sd[j] = sqrt(g->sd[j] / (ngood - 1));
      /*printf("\r");*/
   }
   /*printf("\n");*/

   free(tmp);

   return SUCCESS;
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
   int i, n = 0, p = 0, len, k;
   char *filename_bin = NULL,
        *filename_scale = "scale.bin",
	*filename_beta_out = NULL,
	*filename_folds_ind = NULL;
   short encoded = TRUE,
	 binformat = BINFORMAT_PLINK;
   gmatrix g;
   char tmp[100];

   for(i = 1 ; i < argc ; i++)
   {
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
      else if(strcmp2(argv[i], "-foldind"))
      {
	 i++;
	 filename_folds_ind = argv[i];
      }
   }

   if(!filename_bin || n == 0 || p == 0)
   {
      printf("scale: -bin <filein> [-scale <fileout>] \
[-betafile <betafile>] -n #n -p #p \
[-foldind <folds ind file>]\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_bin, n, p,
	    NULL, YFORMAT01, MODEL_LINEAR, MODELTYPE_REGRESSION,
	    encoded, binformat, filename_folds_ind,
	    MODE_TRAIN, NULL, NULL, NULL))
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
     
   if(filename_beta_out)
      free(filename_beta_out);

   return EXIT_SUCCESS;
}

