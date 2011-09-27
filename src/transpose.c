#include <stdio.h>
#include <stdlib.h>
#include "cd.h"

int transpose(char *filename_in, char *filename_out, const int n,
   const int p, const unsigned int bufsize)
{
   long i, j = 0, bufs = 0;
   unsigned long k, arrlen;
   FILE *in = NULL,
        *out = NULL;
   dtype **buf = NULL,
         *arr = NULL;
   size_t *ret;

   MALLOCTEST(buf, sizeof(dtype*) * n)
   for(i = 0 ; i < n ; i++)
      MALLOCTEST(buf[i], sizeof(dtype) * bufsize)

   MALLOCTEST(ret, sizeof(size_t) * n)
   MALLOCTEST(arr, sizeof(dtype) * bufsize * n);

   FOPENTEST(in, filename_in, "rb");
   FOPENTEST(out, filename_out, "wb");

   while(j < p)
   {
      /* read an (n by bufsize) block */
      for(i = 0 ; i < n ; i++)
      {
	 /* seek to first variable in the block */
	 FSEEKOTEST(in, (unsigned long long)(p * i + bufsize * bufs), SEEK_SET)
	 ret[i] = fread(buf[i], sizeof(dtype),
	       fminl(bufsize, p - bufsize * bufs), in);

	 if(ret[i] == 0 && !feof(in))
	 {
	    fprintf(stderr, "error %d reading file %s",
		  ferror(in), filename_in);
	    return FAILURE;
	 }
      }

      /* unroll block into contiguous array */
      arrlen = 0;
      for(k = 0 ; k < bufsize ; k++)
      {
	 for(i = 0 ; i < n ; i++)
	    if(k < ret[i])
	    {
	       arr[arrlen] = buf[i][k];
	       arrlen++;
	    }
      }

      FWRITETEST(arr, sizeof(dtype), arrlen, out);

      j += bufsize;
      bufs++;
      printf("%lu\n", j);
   }
   printf("\n");

   fclose(in);
   fflush(out);
   fclose(out);
   for(i = 0 ; i < n ; i++)
      free(buf[i]);
   free(buf);
   free(ret);
   free(arr);

   return SUCCESS;
}

/* Transpose a row-wise file to a column-wise file 
 *
 * transpose is compression/coding agnostic, as long as the original matrix
 * row structure is kept intact.
 *
 * The number of columns p should be *all* columns, i.e., inclusive of
 * the y column
 * */
int main(int argc, char *argv[])
{
   int i, n = 0, p = 0;
   char *filename_in = NULL,
        *filename_out = NULL;

   /* bufsize is measured in number of variables, so it's size
    * is sizeof(dtype) * bufsize * n bytes*/
   int bufsize = 0;

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
      else if(strcmp2(argv[i], "-bufsize"))
      {
	 i++;
	 bufsize = (int)atol(argv[i]);
      }
   }

   if(!filename_in || !filename_out || n == 0 || p == 0)
   {
      printf("usage: transpose -fin <filename_in> -fout <filename_out> \
-n #n -p #p\n");
      return EXIT_FAILURE;
   }

   /* Will use about 256MB of RAM: 128MB of buffer + 128MB for output array */
   if(bufsize == 0)
      bufsize = fminl(134217728 * sizeof(dtype) / n, p);

   /* don't forget y is a row too */
   if(!transpose(filename_in, filename_out, n, p, bufsize))
      return EXIT_FAILURE;

   free(filename_in);
   free(filename_out);

   return EXIT_SUCCESS;
}

