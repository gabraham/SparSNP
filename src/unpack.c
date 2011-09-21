/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include "common.h"
#include "gmatrix.h"

int unpack(gmatrix *g, char *filename_out)
{
   int i, j, p1;
   sample sm;
   FILE *out = NULL;
   void *tmp = NULL;
   int start;
   int size = g->scalefile ? sizeof(double) : sizeof(char);

   if(!sample_init(&sm))
      return FAILURE;

   MALLOCTEST(tmp, size * g->n);
   FOPENTEST(out, filename_out, "wb");

   printf("unpack: size=%d scalefile=%s\n", size, g->scalefile);

   if(g->binformat == BINFORMAT_BIN)
   {
      if(g->scalefile)
      {
	 for(i = g->n - 1; i >= 0 ; --i)
	    *((double*)tmp + i) = (double)g->y[i];
      }
      else
      {
	 for(i = g->n - 1; i >= 0 ; --i)
	    *((char*)tmp + i) = (char)g->y[i];
      }
      FWRITETEST(tmp, size, g->n, out);
      start = 1;
      p1 = g->p + 1;
   }
   else
   {
      start = 1;
      p1 = g->p + 1;
   }

   /* ignore intercept */
   for(j = start ; j < p1 ; j++)
   {
      g->nextcol(g, &sm, j, NA_ACTION_ZERO);
      if(g->scalefile)
      {
	 for(i = g->n - 1 ; i >= 0 ; --i)
	    *((double*)tmp + i) = (double)sm.x[i];
      }
      else
      {
	 for(i = g->n - 1 ; i >= 0 ; --i)
	    *((char*)tmp + i) = (char)sm.x[i];
      }
      FWRITETEST(tmp, size, g->n, out);
   }

   fflush(out);
   fclose(out);

   free(tmp);

   return SUCCESS;
}

int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *filename_bin = NULL,
	*filename_out = NULL,
	*file_scale = NULL;
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
      else if(strcmp2(argv[i], "-plink"))
      {
	 binformat = BINFORMAT_PLINK;
      }
      else if(strcmp2(argv[i], "-scale"))
      {
	 i++;
	 file_scale = argv[i];
      }
   }

   if(!filename_bin || !filename_out || n == 0 || p == 0)
   {
      printf("usage: unpack -bin <binfile> -out <outfile> \
-n <#n> -p <#p>\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_bin, n, p, file_scale,
	 YFORMAT01, MODEL_LINEAR, MODELTYPE_REGRESSION, TRUE, binformat,
	 NULL, MODE_TRAIN, NULL, NULL, NULL))
      return EXIT_FAILURE;

   if(!unpack(&g, filename_out))
   {
      gmatrix_free(&g);
      return EXIT_FAILURE;
   }

   gmatrix_free(&g);
   return EXIT_SUCCESS;
}

