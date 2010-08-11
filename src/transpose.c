#include <stdio.h>
#include <stdlib.h>
#include "cd.h"

int transpose(char *filename_in, char *filename_out, int n, int p)
{
   int i, j;
   FILE *in = NULL, *out = NULL;
   dtype *buf = NULL;

   MALLOCTEST(buf, sizeof(dtype) * n);

   FOPENTEST(in, filename_in, "rb");
   FOPENTEST(out, filename_in, "wb");

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	/* FREADTEST(buf[i]) */
      }
   }

   fclose(in);
   fclose(out)
   free(buf);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file */
int main(int argc, char *argv[])
{
   int i, n, p;
   char *filename_in, filename_out;

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
   }

   /* don't forget y is a row too */
   if(!transpose(filename_in, filename_out, n, p + 1))
      return EXIT_FAILURE;

   return EXIT_SUCCESS;
}

