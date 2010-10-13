#include <stdio.h>
#include <stdlib.h>
#include "cd.h"

/* Converts hapgen output to plink PED
 *
 * Assumes:
 * - hapgen data is in ASCII, genotype encoded as one of 0/1/2
 * - after each genotype there is a space
 * - after the final space there is a new line
 */
int hapgen2ped(char *x_filename_in, char *y_filename_in, char *filename_out,
   const int n, const int p, const int bufsize)
{
   int i, j, k, m, rem, bufctr;
   FILE *x_in = NULL, *y_in = NULL, *out = NULL;
   char **buf_in = NULL,
        **buf_out = NULL,
        **buf_y = NULL;
   const int p2 = 2 * p + 1;
   const short ASCII_DIFF = 48;
   char const *restrict ptr = NULL;

   /* The spaces are important */
   const char *alleles[3] = {
      "A A ", /* 0 */
      "A B ", /* 1 */
      "B B "  /* 2 */
   };
   const char constfields[4] = {
      '1', /* Family ID */
      '0', /* Paternal ID */
      '0', /* Maternal ID */
      '1'  /* Sex */
   };
   const char affected[2] = {
      '1', /* unaffected */
      '2'  /* affected */
   };


   /* make room for spaces as well, including the EOL*/
   MALLOCTEST(buf_in, sizeof(char*) * bufsize)
   MALLOCTEST(buf_out, sizeof(char*) * bufsize)
   MALLOCTEST(buf_y, sizeof(char*) * bufsize)

   for(k = 0 ; k < bufsize ; k++)
   {
      MALLOCTEST(buf_in[k], sizeof(char) * p2);
      MALLOCTEST(buf_out[k], sizeof(char) * (4 * p + 1));
      MALLOCTEST(buf_y[k], sizeof(char) * 2);
      buf_out[k][4 * p] = '\0';
   }

   FOPENTEST(x_in, x_filename_in, "rt");
   FOPENTEST(y_in, y_filename_in, "rt");
   FOPENTEST(out, filename_out, "wt");

   /* read data in blocks */
   bufctr = 0;
   for(i = 0 ; i < n ; )
   {
      printf("%d of %d", i, n);
      fflush(stdout);

      rem = FMIN(bufsize, n - bufctr * bufsize);
      for(k = 0 ; k < rem ; k++)
      {
	 FREADTEST(buf_y[k], sizeof(char), 2, y_in);
      	 FREADTEST(buf_in[k], sizeof(char), p2, x_in);

      	 for(j = 0 ; j < p ; j++)
      	 {
      	    ptr = alleles[buf_in[k][2 * j] - ASCII_DIFF];
	    m = 4 * j;
      	    buf_out[k][m]     = ptr[0];
      	    buf_out[k][m + 1] = ptr[1];
      	    buf_out[k][m + 2] = ptr[2];
      	    buf_out[k][m + 3] = ptr[3];
      	 }
      }
	 
      for(k = 0 ; k < rem ; k++)
      {
	 fprintf(out, "%c %d %c %c %c %c ", constfields[0], i + k + 1,
	       constfields[1], constfields[2], constfields[3],
	       affected[buf_y[k][0] - ASCII_DIFF]);
      	 fprintf(out, "%s\n", buf_out[k]);
      }

      i += bufsize;
      bufctr++;
      printf("\r");
   }
   printf("\n");

   fclose(x_in);
   fclose(y_in);
   fflush(out);
   fclose(out);
   for(k = 0 ; k < bufsize ; k++)
   {
      free(buf_in[k]);
      free(buf_out[k]);
      free(buf_y[k]);
   }
   free(buf_in);
   free(buf_out);
   free(buf_y);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file */
int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *x_filename_in = NULL,
        *y_filename_in = NULL,
	*filename_out = NULL;

   int bufsize = 0;

   MALLOCTEST(x_filename_in, sizeof(char) * MAX_STR_LEN)
   MALLOCTEST(y_filename_in, sizeof(char) * MAX_STR_LEN)
   MALLOCTEST(filename_out, sizeof(char) * MAX_STR_LEN)

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-finx"))
      {
	 i++;
	 strncpy(x_filename_in, argv[i], MAX_STR_LEN);
      }
      if(strcmp2(argv[i], "-finy"))
      {
	 i++;
	 strncpy(y_filename_in, argv[i], MAX_STR_LEN);
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
      else if(strcmp2(argv[i], "-bufsize"))
      {
	 i++;
	 bufsize = (int)atol(argv[i]);
      }
   }

   if(!x_filename_in || !y_filename_in || !filename_out || n == 0 || p == 0)
   {
      printf("usage: hapgen2ped -finx <x_filename_in> -finy <y_filename_in> \
-fout <filename_out> -n #n -p #p\n");
      return EXIT_FAILURE;
   }

   /* bufsize is a multiple of rows (2 * p + 1 values)  */
   if(bufsize == 0)
      bufsize = fminl(3, n);

   /* don't forget y is a row too */
   if(!hapgen2ped(x_filename_in, y_filename_in,
	 filename_out, n, p, bufsize))
      return EXIT_FAILURE;

   free(x_filename_in);
   free(y_filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

