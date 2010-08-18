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
   const unsigned int n, const unsigned int p, const unsigned int bufsize)
{
   unsigned long i, j;
   FILE *x_in = NULL, *y_in = NULL, *out = NULL;
   char *buf_in = NULL, *buf_out = NULL;
   char y[2];
   const unsigned int p1 = p + 1, p2 = 2 * p + 1;
   const short ASCII_DIFF = 48;
   char *ptr = NULL;

   /* The spaces are important */
   char *alleles[3] = {
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
   MALLOCTEST(buf_in, sizeof(char) * p2)
   MALLOCTEST(buf_out, sizeof(char) * (4 * p + 1))

   FOPENTEST(x_in, x_filename_in, "rt");
   FOPENTEST(y_in, y_filename_in, "rt");
   FOPENTEST(out, filename_out, "wt");

   buf_out[4 * p] = '\0';

   for(i = 0 ; i < n ; i++)
   {
      printf("%ld\r", i);
      FREADTEST(&y, sizeof(char), 2, y_in)
      
      fprintf(out, "%c %lu %c %c %c %c ", constfields[0], i + 1, constfields[1],
	    constfields[2], constfields[3], affected[y[0] - ASCII_DIFF]);

      FREADTEST(buf_in, sizeof(char), p2, x_in)

      for(j = 0 ; j < p ; j++)
      {
	 ptr = alleles[buf_in[2 * j] - ASCII_DIFF];
	 buf_out[4 * j] = ptr[0];
	 buf_out[4 * j + 1] = ptr[1];
	 buf_out[4 * j + 2] = ptr[2];
	 buf_out[4 * j + 3] = ptr[3];
      }
	 
      fprintf(out, "%s\n", buf_out);
   }
   printf("\n");

   fclose(x_in);
   fclose(y_in);
   fflush(out);
   fclose(out);
   free(buf_in);
   free(buf_out);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file */
int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *x_filename_in = NULL,
        *y_filename_in = NULL,
	*filename_out = NULL;

   unsigned int bufsize = 0;

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

   /* bufsize is a multiple of rows (p+1 values)  */
   if(bufsize == 0)
      bufsize = fminl(134217728 * sizeof(dtype) / (p + 1), n);

   /* don't forget y is a row too */
   if(!hapgen2ped(x_filename_in, y_filename_in,
	 filename_out, n, p, bufsize))
      return EXIT_FAILURE;

   free(x_filename_in);
   free(y_filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

