#include  <time.h>
#include "common.h"
#include "gmatrix.h"

int subsample(gmatrix *g, char *filename_out, int *subvec, int nsub)
{
   int i, k, j, p1 = g->p + 1;
   const int numencb = (int)ceil(nsub / (double)PACK_DENSITY);
   sample sm;
   FILE *out;
   unsigned char *tmp = NULL;
   unsigned char *enc_buf = NULL;

   if(!sample_init(&sm))
      return FAILURE;

   MALLOCTEST(tmp, sizeof(char) * nsub);
   MALLOCTEST(enc_buf, sizeof(unsigned char) * numencb);
   FOPENTEST(out, filename_out, "wb");

   printf("Using %d samples, %d variables (exc. intercept)\n", nsub, g->p);
   k = 0;
   for(i = 0 ; i < g->n ; i++)
      if(subvec[i])
	 tmp[k++] = (unsigned char)g->y[i];

   encode(enc_buf, tmp, nsub);
   FWRITETEST(enc_buf, sizeof(unsigned char), numencb, out);
   
   /* read a variable and write it, ignore intercept */
   for(j = 1 ; j < p1 ; j++)
   {
      g->nextcol(g, &sm, j, NA_ACTION_RANDOM);

      k = 0;
      for(i = 0 ; i < g->n ; i++)
	 if(subvec[i])
	    tmp[k++] = (unsigned char)sm.x[i];
      
      encode(enc_buf, tmp, nsub);
      FWRITETEST(enc_buf, sizeof(unsigned char), numencb, out);
   }

   fflush(out);
   fclose(out);

   free(tmp);
   free(enc_buf);

   return SUCCESS;
}

/* Populate the array subvec with random 0/1, so that the
 * total number of 1s is roughly ns (out of n).
 */
int init_subsamples(int *subvec, int n, int ns, double *y)
{
   int i;
   int num1 = 0; /* number of 1s in subvec */
   long seed = time(NULL);
   double prop = (double)ns / n;
   int monoclass = TRUE;
   int ny1 = 0;

   srand48(seed);

   while(monoclass)
   {
      num1 = 0;
      ny1 = 0;
      for(i = 0 ; i < n ; i++)
      {
         subvec[i] = (drand48() < prop);
         num1 += subvec[i];
	 ny1 += subvec[i] && y[i];
      }

      monoclass = (ny1 == 0 || ny1 == n);
   }

   printf("Class balance (0/1): %d/%d, total=%d\n", num1-ny1, ny1, num1);

   return num1;
}

int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *filename_bin = NULL,
	*filename_out = NULL;
   short binformat = BINFORMAT_BIN;
   gmatrix g;
   int *subvec = NULL;
   int ns = 0,
       nsemp = 0;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-bin"))
      {
	 i++;
	 filename_bin = argv[i];
      }
      if(strcmp2(argv[i], "-out"))
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
      else if(strcmp2(argv[i], "-ns"))
      {
	 i++;
	 ns = (int)atoi(argv[i]);
      }
   }

   if(!filename_bin || !filename_out || n == 0 || p == 0)
   {
      printf("usage: subsample -bin <binfile> -out <outfile> \
-n <#n> -p <#p> -ns <#num subsample>\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_bin, n, p, NULL,
	 YFORMAT01, MODEL_LINEAR, MODELTYPE_REGRESSION, TRUE, binformat,
	 NULL, MODE_TRAIN, NULL, NULL, NULL))
      return EXIT_FAILURE;

   MALLOCTEST2(subvec, sizeof(int) * n);

   /* clip to reasonable range */
   ns = (ns <= 1) ? 1 : (ns >= n - 1 ? n - 1 : ns);
   nsemp = init_subsamples(subvec, n, ns, g.y);

   if(!subsample(&g, filename_out, subvec, nsemp))
   {
      gmatrix_free(&g);
      return EXIT_FAILURE;
   }

   gmatrix_free(&g);
   FREENULL(subvec);

   return EXIT_SUCCESS;
}

