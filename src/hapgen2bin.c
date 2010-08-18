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
   unsigned long i, j = 0, k, bufs = 0, arrlen = 0, len;
   FILE *x_in = NULL,
        *y_in = NULL,
	*out = NULL;
   unsigned char *buf_hg = NULL,
		 **buf_bin = NULL,
		 **buf_bin_enc = NULL,
		 *arr = NULL,
		 *buf = NULL,
		 **buf_ptr = NULL,
		 *y_bin = NULL,
		 *y_bin_enc = NULL;
   const unsigned int p1 = p + 1, /* genotypes without spaces, and EOL */
		      p2 = 2 * p + 1, /* genotypes with spaces, and EOL */
		      ncoded = (unsigned int)
			    (ceil(n / (double)PACK_DENSITY));
   size_t *nread = NULL;
   const unsigned char ASCII_DIFF = 48;

   MALLOCTEST(buf_hg, sizeof(unsigned char) * bufsize);

   /* n buffers of size bufsize */
   MALLOCTEST(buf_bin, sizeof(unsigned char*) * n);
   for(i = 0 ; i < n ; i++)
      MALLOCTEST(buf_bin[i], sizeof(unsigned char) * bufsize);

   if(encodeflag)
   {
      /* Encode across the columns. buf_bin_enc is an array of bufsize
       * buffers, each representing one encoded variable of length ncoded
       */
      MALLOCTEST(buf_bin_enc, sizeof(unsigned char*) * bufsize);
      for(i = 0 ; i < ncoded ; i++)
         MALLOCTEST(buf_bin_enc[i], sizeof(unsigned char) * ncoded);

      /* arr is an unrolled buffer, to enable bulk writing to disk */
      MALLOCTEST(arr, sizeof(unsigned char) * bufsize * ncoded);
   }
   else
      MALLOCTEST(arr, sizeof(unsigned char) * bufsize * n);
   
   /* return value from reading x data, used to check we don't try to write
    * too many values to the buffer when EOF happens mid-buffer */
   MALLOCTEST(nread, sizeof(size_t) * n);

   FOPENTEST(x_in, x_filename_in, "rt");
   FOPENTEST(y_in, y_filename_in, "rt");
   FOPENTEST(out, filename_out, "wb");

   /* read y and EOL */
   FREADTEST(buf_hg, sizeof(unsigned char), bufsize, y_in);
   for(i = 0 ; i < n ; i++)
      y_bin[i] = buf_hg[2 * i] - ASCII_DIFF;

   /* process y */
   if(encodeflag)
   {
      MALLOCTEST(y_bin_enc, sizeof(unsigned char) * ncoded);
      encode(y_bin_enc, y_bin, n);
      FWRITETEST(y_bin_enc, sizeof(unsigned char), ncoded, out);
      free(y_bin_enc);
      y_bin_enc = NULL;
   }
   else
      FWRITETEST(y_bin, sizeof(unsigned char), n, out);

   free(y_bin);
   y_bin = NULL;

   /* process the x variables */
   while(j < p)
   {
      for(i = 0 ; i < n ; i++)
      {
         printf("%ld\r", i);
 
	 /* seek to first variable in the (n by bufsize) buffer, then fill the
	  * buffer as much as possible */
	 FSEEKOTEST(x_in, (unsigned long long)(p2 * i + bufsize * bufs),
	       SEEK_SET);
	 nread[i] = fread(buf_hg, sizeof(unsigned char),
	       fminl(bufsize, p2 - bufsize * bufs), x_in);
 
	 if(nread[i] == 0 && !feof(x_in))
	 {
	    fprintf(stderr, "error %d reading file %s",
		  ferror(x_in), x_filename_in);
	    return FAILURE;
	 }

	 /* convert to genotypes, dropping spaces */
         for(k = 0 ; k < p ; k++)
	    buf_bin[i][k + 1] = (dtype)(buf_hg[2 * k] - ASCII_DIFF);
   
	 arrlen = 0;
	 buf_ptr = buf_bin;
	 len = n;
	 
         /* encode each column *separately*, slighly space-inefficient
          * but easier to decode */ 
	 if(encodeflag)
	 {
	    for(k = 0 ; k < bufsize ; k++)
	       encode(buf_bin_enc[k], buf_bin[k], p1);

	    buf_ptr = buf_bin_enc;
	    len = ncoded;
	 }
  
	 /* unroll block into contiguous array */
	 for(k = 0 ; k < bufsize ; k++)
	    for(i = 0 ; i < len ; i++)
	       if(k < nread[i])
	       {
		  arr[arrlen] = buf_ptr[i][k];
		  arrlen++;
	       }

	 FWRITETEST(arr, sizeof(dtype), arrlen, out);

	 j += bufsize;
	 bufs++;

      }
   }
   printf("\n");

   fclose(x_in);
   fclose(y_in);
   fflush(out);
   fclose(out);
   free(buf_hg);

   for(i = 0 ; i < n ; i++)
      free(buf_bin[i]);
   free(buf_bin);

   if(buf_bin_enc)
   {
      for(i = 0 ; i < ncoded ; i++)
         free(buf_bin_enc[i]);
      free(buf_bin_enc);
   }

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

   /* bufsize is a multiple of rows (p+1 values) */
   if(bufsize == 0)
      bufsize = fminl(134217728 * sizeof(dtype) / p, n);

   /* don't forget y is a row too */
   if(!hapgen2bin(x_filename_in, y_filename_in,
	 filename_out, n, p, bufsize, encodeflag))
      return EXIT_FAILURE;

   free(x_filename_in);
   free(y_filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

