#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common.h"

/* Concatenates *text* files by columns 
 */
int cbind(const int nfiles, char **filenames_in,
      const char *filename_out, const int n,
      const int p)
{
   int i, j, p2 = 2 * p, p22 = p2 + 2;
   FILE **in = NULL, *out = NULL;
   char *buf = NULL;

   MALLOCTEST(in, sizeof(FILE) * nfiles);
   MALLOCTEST(buf, sizeof(char) * p22);
   FOPENTEST(out, filename_out, "wt");

   for(j = 0 ; j < nfiles ; j++)
      FOPENTEST(in[j], filenames_in[j], "rt");

   for(i = 0 ; i < n ; i++)
   {
      printf("%d of %d", i, n);
      for(j = 0 ; j < nfiles ; j++)
      {
	 if(!(buf = fgets(buf, p22, in[j])))
	 {
	    fprintf(stderr, "error reading from input file");
	    return FAILURE;
	 }
	 buf[strlen(buf) - 1] = '\0';
	 fprintf(out, "%s", buf);
      }
      fprintf(out, "\n");
      printf("\r");
   }


   for(j = 0 ; j < nfiles ; j++)
      fclose(in[j]);
   free(in);

   fclose(out);
   free(buf);

   return SUCCESS;
}

int main(int argc, char *argv[])
{
   int i, j, k, n = 0, p = 0, nfiles = 0;
   char **filenames_in = NULL,
	 *filename_out = NULL;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-in"))
      {
	 j = ++i;
	 k = 0;
	 while(j < argc && argv[j][0] != '-')
	 {
	    k++;
	    REALLOCTEST(filenames_in, filenames_in, sizeof(char*) * k);
	    MALLOCTEST(filenames_in[k - 1], sizeof(char) * (strlen(argv[j]) + 1));
	    strcpy(filenames_in[k - 1], argv[j]);
	    j++;
	 }
	 if(j < argc && argv[j][0] == '-')
	    i = j - 1;

	 nfiles = k;
      }
      if(strcmp2(argv[i], "-out"))
      {
	 i++;
	 MALLOCTEST(filename_out, sizeof(char) * (strlen(argv[i]) + 1));
	 strcpy(filename_out, argv[i]);
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

   printf("n: %d   p:%d\n", n, p);
   if(!filenames_in || !filename_out || n == 0 || p == 0)
   {
      printf("usage: cbind -in <file_in1> <file_in2> ... -out <file_out> \
-n #n -p #p\n");
      for(i = 0 ; i < nfiles ; i++)
	 free(filenames_in[i]);
      free(filenames_in);
      free(filename_out);
      return EXIT_FAILURE;
   }

   
   if(!(cbind(nfiles, filenames_in, filename_out, n, p)))
   {
      for(i = 0 ; i < nfiles ; i++)
	 free(filenames_in[i]);
      free(filenames_in);
      free(filename_out);
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

