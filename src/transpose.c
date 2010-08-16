#include <stdio.h>
#include <stdlib.h>
#include "cd.h"

int transpose(char *filename_in, char *filename_out, const int n, const int p)
{
   int i, j;
   FILE *in = NULL, *out = NULL;
   dtype *buf = NULL;

   MALLOCTEST(buf, sizeof(dtype) * n);

   FOPENTEST(in, filename_in, "rb");
   FOPENTEST(out, filename_out, "wb");

   for(j = 0 ; j < p ; j++)
   {
      /* read one datum and skip one row of variables */
      for(i = 0 ; i < n ; i++)
      {
	 FSEEKTEST(in, p * i + j, SEEK_SET)
	 FREADTEST(buf + i, sizeof(dtype), 1, in)
      }
      FWRITETEST(buf, sizeof(dtype), n, out)
   }

   fclose(in);
   fflush(out);
   fclose(out);
   free(buf);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file */
int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *filename_in = NULL, *filename_out = NULL;

   MALLOCTEST(filename_in, sizeof(char) * MAX_STR_LEN)
   MALLOCTEST(filename_out, sizeof(char) * MAX_STR_LEN)

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-fin"))
      {
	 i++;
	 strncpy(filename_in, argv[i], MAX_STR_LEN);
      }
      else if(strcmp2(argv[i], "-fout"))
      {
	 i++;
	 strncpy(filename_out, argv[i], MAX_STR_LEN);
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

   if(!filename_in || !filename_out || n == 0 || p == 0)
   {
      printf("usage: transpose -fin <filename_in> -fout <filename_out> \
-n #n -p #p\n");
      return EXIT_FAILURE;
   }



   /* don't forget y is a row too */
   if(!transpose(filename_in, filename_out, n, p + 1))
      return EXIT_FAILURE;

   free(filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

