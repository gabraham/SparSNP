#include "common.h"
#include "gmatrix.h"

int unpack(gmatrix *g, char *filename_out)
{
   int j, p1 = g->p + 1;
   sample sm;
   FILE *out;

   if(!sample_init(&sm, g->n))
      return FAILURE;

   FOPENTEST(out, filename_out, "wb");
   FWRITETEST(g->y, sizeof(double), g->n, out);

   /* ignore intercept */
   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d of %d", j, p1);
      g->nextcol(g, &sm, j);
      FWRITETEST(sm.x, sizeof(double), g->n, out);
      printf("\r");
   }
   printf("\n");

   fflush(out);
   fclose(out);

   return SUCCESS;
}

int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *filename_bin = NULL,
	*filename_out = NULL;
   short binformat = BINFORMAT_BIN;
   gmatrix g;

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
   }

   if(!filename_bin || !filename_out || n == 0 || p == 0)
   {
      printf("usage: unpack -bin <binfile> -out <outfile> \
-n <#n> -p <#p>\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_bin, n, p, NULL,
	 YFORMAT01, MODEL_LINEAR, TRUE, binformat,
	 NULL, MODE_TRAIN, NULL))
      return EXIT_FAILURE;

   if(!unpack(&g, filename_out))
   {
      gmatrix_free(&g);
      return EXIT_FAILURE;
   }

   gmatrix_free(&g);
   return EXIT_SUCCESS;
}

