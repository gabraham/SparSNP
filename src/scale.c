#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"

/*
 * Converts a matrix encoded as chars (one byte per entry), scales it, and
 * saves it as a matrix of doubles (8 bytes per entry).
 *
 * Expects data in column major ordering
 */
int scale2(gmatrix *g, char* filename)
{
   int i, j;
   FILE *fout = NULL, *fin = NULL;
   dtype *tmp;
   double *tmp2;
   double delta;

   MALLOCTEST(g->mean, sizeof(double) * (g->p + 1))
   MALLOCTEST(g->sd, sizeof(double) * (g->p + 1))
   g->mean[0] = 0;
   g->sd[0] = 1;

   MALLOCTEST(tmp, sizeof(dtype) * g->n)
   MALLOCTEST(tmp2, sizeof(double*) * g->n)

   FOPENTEST(fin, g->filename, "rb")

   if(filename)
      FOPENTEST(fout, filename, "w")

   /* read y but do not scale it */
   if(g->encoded) {
      FREADTEST(g->encbuf, sizeof(dtype), g->nencb, fin);
      decode(tmp, g->encbuf, g->nencb);
   } else {
      FREADTEST(tmp, sizeof(dtype), g->n, fin);
   }

   if(filename)
   {
      for(i = 0 ; i < g->n ; i++)
	 tmp2[i] = (double)tmp[i];
      FWRITETEST(tmp2, sizeof(double), g->n, fout);
   }

   /* read the data and scale each variable */
   for(j = 1 ; j < g->p + 1 ; j++)
   {
      printf("%d of %d\r", j, g->p);
      if(g->encoded) {
	 FREADTEST(g->encbuf, sizeof(dtype), g->nencb, fin);
	 decode(tmp, g->encbuf, g->nencb);
      } else {
	 FREADTEST(tmp, sizeof(dtype), g->n, fin);
      }

      g->mean[j] = g->sd[j] = 0;
      for(i = 0 ; i < g->n ; i++)
      {
	 tmp2[i] = (double)tmp[i];
	 delta = tmp2[i] - g->mean[j];
	 g->mean[j] += delta / (i + 1);
	 g->sd[j] += delta * (tmp2[i] - g->mean[j]);
      }

      g->sd[j] = sqrt(g->sd[j] / (g->n - 1));

      if(filename)
      {
	 for(i = 0 ; i < g->n ; i++)
      	 {
      	    tmp2[i] = tmp2[i] - g->mean[j];
      	    if(g->sd[j] > SDTHRESH)
      	       tmp2[i] /= g->sd[j];
      	 }
      	   
      	 FWRITETEST(tmp2, sizeof(double), g->n, fout);
      }
   }
   printf("\n");

   fclose(fin);
   if(filename)
      fclose(fout);

   free(tmp);
   free(tmp2);
      
   return SUCCESS;
}

int scale(gmatrix *g, char* filename)
{
   unsigned int i, j, p1 = g->p + 1;
   sample sm;
   double *tmp = NULL, delta;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   MALLOCTEST(g->mean, sizeof(double) * p1);
   MALLOCTEST(g->sd, sizeof(double) * p1);
   g->mean[0] = 0;
   g->sd[0] = 1;

   MALLOCTEST(tmp, sizeof(double) * g->n);
   
   /* read intercept and ignore it*/
   g->nextcol(g, &sm);

   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d of %d", j, p1);
      g->nextcol(g, &sm);
      g->mean[j] = g->sd[j] = 0;
      for(i = 0 ; i < g->n ; i++)
      {
	 delta = sm.x[i] - g->mean[j];
	 g->mean[j] += delta / (i + 1);
	 g->sd[j] += delta * (sm.x[i] - g->mean[j]);
      }

      g->sd[j] = sqrt(g->sd[j] / (g->n - 1));
      printf("\r");
   }
   printf("\n");

   free(tmp);
   sample_free(&sm);

   return SUCCESS;
}

int writescale(char* filename, double *mean, double *sd,
      int p, short ascii)
{
   FILE *out;

   if(ascii)
   {
      if(!writevectorf("mean.csv", mean, p))
	 return FAILURE;

      if(!writevectorf("sd.csv", sd, p))
	 return FAILURE;

      return SUCCESS;
   }
   
   FOPENTEST(out, filename, "wb")
   FWRITETEST(mean, sizeof(double), p, out);
   FWRITETEST(sd, sizeof(double), p, out);
   fclose(out);
   return SUCCESS;
}

int main(int argc, char* argv[])
{
   int i, n = 0, p = 0;
   char *filename_in = NULL,
        *filename_out = NULL,
        *filename_scale = "scale.bin";
   short inmemory = FALSE,
	 ascii = FALSE,
	 encoded = FALSE,
	 binformat = BINFORMAT_BIN;
   gmatrix g;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-fin"))
      {
	 i++;
	 filename_in = argv[i];
      }
      else if(strcmp2(argv[i], "-fout"))
      {
	 i++;
	 filename_out = argv[i];
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
      else if(strcmp2(argv[i], "-inmemory"))
	 inmemory = TRUE;
      else if(strcmp2(argv[i], "-ascii"))
	 ascii = TRUE;
      else if(strcmp2(argv[i], "-encoded"))
	 encoded = TRUE;
      else if(strcmp2(argv[i], "-plink"))
	 binformat = BINFORMAT_PLINK;
   }

   if(filename_in == NULL || n == 0 || p == 0)
   {
      printf("scale: -fin <filein> [-fout <fileout>] \
[-fscale <scalefile>] [-encoded] [-plink] -n #n -p #p\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_in, n, p,
	 inmemory, NULL, YFORMAT01, MODE_TRAIN, encoded, binformat))
      return EXIT_FAILURE;

   if(!scale(&g, filename_out))
      return EXIT_FAILURE;

   if(!writescale(filename_scale, g.mean, g.sd, p + 1, ascii))
      return EXIT_FAILURE;
  
   gmatrix_free(&g);

   return EXIT_SUCCESS;
}

