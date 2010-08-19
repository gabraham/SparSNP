#include <stdio.h>
#include <stdlib.h>
#include "cd.h"
#include "coder.h"

/* Assumes:
 * - hapgen data is in ASCII, genotype encoded as one of 0/1/2/NA
 * - after each genotype there is a space
 * - after the final space there is a new line
 * - hapgen phenotype data is one character per line plus a newline
 * - byte coding: encodes the rows (p+1 bytes)
 */
int hapgen2bin(char *x_filename_in, char *y_filename_in, char *filename_out,
   const unsigned int n, const unsigned int p, const unsigned int bufsize,
   const short encodeflag)
{
   const unsigned int ASCII_SHIFT = 48;
   unsigned int i, j, k, nbufs, nbufb,
	        numencb = (unsigned int)ceil(n / (double)PACK_DENSITY),
		p21 = 2 * p + 1;
   FILE *in_x = NULL,
	*in_y = NULL,
	*out = NULL;
   unsigned char *y_buf = NULL,
	         *y_bin = NULL,
		 *buf = NULL,
		 *enc_buf = NULL,
		 *hg_buf = NULL;
   unsigned long long skip;

   FOPENTEST(in_x, x_filename_in, "rt");
   FOPENTEST(in_y, y_filename_in, "rt");
   FOPENTEST(out, filename_out, "wb");

   MALLOCTEST(y_buf, sizeof(unsigned char) * 2 * n);
   MALLOCTEST(y_bin, sizeof(unsigned char) * n);
   MALLOCTEST(hg_buf, sizeof(unsigned char*) * n * bufsize);
   MALLOCTEST(enc_buf, sizeof(unsigned char) * numencb * bufsize);
   MALLOCTEST(buf, sizeof(unsigned char) * bufsize * 2);

   /* process y */
   FREADTEST(y_buf, sizeof(unsigned char), 2 * n, in_y);
   for(i = 0 ; i < n ; i++)
      y_bin[i] = y_buf[2 * i] - ASCII_SHIFT;

   if(encodeflag) {
      encode(enc_buf, y_bin, n);
      FWRITETEST(enc_buf, sizeof(unsigned char), numencb, out);
   } else {
      FWRITETEST(y_bin, sizeof(unsigned char), n, out);
   }

   free(y_bin);
   free(y_buf);
   y_bin = y_buf = NULL;

   nbufs = 0;
   for(j = 0 ; j < p ; j += bufsize)
   {
      fflush(stdout);
      nbufb = (unsigned int)fmin(bufsize, p - bufsize * nbufs);
      for(i = 0 ; i < n ; i++)
      {
         /* read one block of genotypes, skipping over spaces and EOL */
	 skip = (unsigned long long)i * p21 
	       + (unsigned long long)j * 2;
         FSEEKOTEST(in_x, skip, SEEK_SET);
         FREADTEST(buf, sizeof(unsigned char), 2 * nbufb, in_x);
	 fflush(stdout);
	 for(k = 0 ; k < nbufb ; k++)
	    hg_buf[k * n + i] = buf[2 * k] - ASCII_SHIFT;
      }

      /* encode */
      if(encodeflag) {
	 /* encode each variable separately to prevent two falling into the
	  * same encoded byte */
	 for(k = 0 ; k < nbufb ; k++)
	    encode(enc_buf + k * numencb, hg_buf + k * n, n);
	 FWRITETEST(enc_buf, sizeof(unsigned char), numencb * nbufb, out);
      } else {
	 FWRITETEST(hg_buf, sizeof(unsigned char), n * nbufb, out);
      }
      nbufs++;
   }

   fclose(in_x);
   fclose(in_y);
   fclose(out);
   free(enc_buf);
   free(hg_buf);
   free(buf);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file */
int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *x_filename_in = NULL,
        *y_filename_in = NULL,
	*filename_out = NULL;

   short encodeflag = FALSE;

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
      else if(strcmp2(argv[i], "-encode"))
      {
	 encodeflag = TRUE;
      }
   }

   if(!x_filename_in || !y_filename_in || !filename_out || n == 0 || p == 0)
   {
      printf("usage: hapgen2bin -finx <x_filename_in> -finy <y_filename_in> \
-fout <filename_out> -n #n -p #p [-encode]\n");
      return EXIT_FAILURE;
   }

   if(bufsize == 0)
      bufsize = 1024;
   else if(bufsize > p)
      bufsize = p;

   /* don't forget y is a row too */
   if(!hapgen2bin(x_filename_in, y_filename_in,
	 filename_out, n, p, bufsize, encodeflag))
      return EXIT_FAILURE;

   free(x_filename_in);
   free(y_filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

