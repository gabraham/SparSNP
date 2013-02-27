/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include "common.h"
#include "sparsnp.h"

#define UNPACK_FORMAT_CHAR 1
#define UNPACK_FORMAT_DOUBLE 2

int unpack(gmatrix *g, char *filename_out, int dim, int format)
{
   int i, j, p1;
   sample sm;
   FILE *out = NULL;
   void *tmp = NULL;
   int start;
   /*int size = g->scalefile ? sizeof(double) : sizeof(char);*/
   int size = (format == UNPACK_FORMAT_CHAR) ? sizeof(char) : sizeof(double);

   if(!sample_init(&sm))
      return FAILURE;

   MALLOCTEST(tmp, size * g->n);
   FOPENTEST(out, filename_out, "wb");

   printf("unpack: size=%d scalefile=%s\n", size, g->scalefile);
   start = 1;
   p1 = g->p + 1;

   /* Write dimensions of the matrix */
   if(dim)
   {
      printf("unpack: n=%d p=%d\n", g->n, g->p);
      FWRITETEST(&(g->n), sizeof(g->n), 1, out);
      FWRITETEST(&(g->p), sizeof(g->p), 1, out);
   }

   /* ignore intercept, start from j=1 */
   for(j = start ; j < p1 ; j++)
   {
      if(!g->nextcol(g, &sm, j, NA_ACTION_PROPORTIONAL))
	 return FAILURE;

      if(format == UNPACK_FORMAT_DOUBLE)
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
   int dim = TRUE;
   char *filename_bed = NULL;
   char *filename_out = NULL;
   char *file_scale = NULL;
   gmatrix g;
   int unpack_format = UNPACK_FORMAT_DOUBLE;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-bed"))
      {
	 i++;
	 filename_bed = argv[i];
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
      else if(strcmp2(argv[i], "-scale"))
      {
	 i++;
	 file_scale = argv[i];
      }
      else if(strcmp2(argv[i], "-dim"))
      {
	 i++;
	 dim = TRUE;
      }
      else if(strcmp2(argv[i], "-format"))
      {
	 i++;
	 unpack_format = (int)atof(argv[i]);
      }
   }

   if(!filename_bed || !filename_out || n == 0 || p == 0)
   {
      printf("usage: unpack -bed <bedfile> -out <outfile> \
-n <#n> -p <#p> [-format double|char] [-dim]\n");
      return EXIT_FAILURE;
   }

   printf("unpacking %s to file %s\n", filename_bed, filename_out);

   if(!gmatrix_init(&g, filename_bed, n, p, file_scale,
	 YFORMAT01, 0, MODEL_LINEAR, MODELTYPE_REGRESSION, TRUE,
	 NULL, MODE_TRAIN, NULL, NULL, FALSE, FALSE, 0, 0, FALSE,
	 CACHE_MEM_DEFAULT))
      return EXIT_FAILURE;

   if(!unpack(&g, filename_out, dim, unpack_format))
   {
      gmatrix_free(&g);
      return EXIT_FAILURE;
   }

   gmatrix_free(&g);
   return EXIT_SUCCESS;
}

