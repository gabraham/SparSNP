#include <stdlib.h>
#include <stdio.h>
#include "coder.h"

/*
 * Converts an array of short integers {0,1,2,3} (3 is NA) to binary 
 * {00, 01, 10, 11}, so that each byte contains 4 values
 *
 * length of out should be ceiling(length(in) / 4.0), so it should be padded
 * with zeros if the number of bytes is not a multiple of 4. This function
 * does not check for alignment!
 * 
 * n: length of in
 *
 */
void encode(unsigned char *out, const unsigned char *in, const unsigned int n)
{
   unsigned int i, j, k = 0;

   for(i = 0 ; i < n ; i += PACK_DENSITY)
   {  
      printf("n:%d i:%d k:%d\n", n, i, k);
      out[k] = 0;
      for(j = 0 ; j < PACK_DENSITY ; j++)
	 if(i + j < n)
	    out[k] |= in[i + j] << (2 * j);
      k++;
   }
}

/*
 * n: length of in
 */
void decode(unsigned char *out, const unsigned char *in, const unsigned int n)
{
   unsigned int i, j;
   unsigned char val, tmp;
   const unsigned char masks[PACK_DENSITY] = {
      3 << 2 * 0, /* 3 is 11 in binary */
      3 << 2 * 1,
      3 << 2 * 2,
      3 << 2 * 3,
   };

   for(i = 0 ; i < n ; i++)
   {
      tmp = in[i];
      for(j = 0; j < PACK_DENSITY ; j++)
      {
	 val = (in[i] & masks[j]) >> (2 * j);
	 out[PACK_DENSITY * i + j] = val; 
	 tmp -= val;
      }
   }
}

/* int main(void)
{
   int i;
   const unsigned char in[12] = {
      0, 1, 2, 3,
      1, 0, 2, 0,
      1, 3, 0, 1
   };
   unsigned char out[3] = {0, 0, 0};
   unsigned char out2[12] = {
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0
   };

   encode(out, in, 12);

   for(i = 0 ; i < 3 ; i++)
      printf("encode: %d\n", out[i]);
   printf("\n");

   decode(out2, out, 3);

   for(i = 0 ; i < 12 ; i++)
      printf("decode: %d (%d)\n", out2[i], in[i]);
   printf("\n");


   
   return EXIT_SUCCESS;
}
*/

