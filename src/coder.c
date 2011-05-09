#include <stdlib.h>
#include <stdio.h>
#include "coder.h"

/*
 * Converts an array of short integers {0,1,2,3} (3 is NA) to binary 
 * {00, 01, 10, 11}, so that each byte contains 4 genotypes.
 *
 * length of out should be ceiling(length(in) / 4.0), so it should be padded
 * with zeros if the number of bytes is not a multiple of 4. This function
 * does not check for alignment!
 * 
 * n: length of "in" buffer
 *
 */
void encode(unsigned char *out, const unsigned char *in, const int n)
{
   int i, j, k = 0;

   for(i = 0 ; i < n ; i += PACK_DENSITY)
   {  
      out[k] = 0;
      for(j = 0 ; j < PACK_DENSITY ; j++)
	 if(i + j < n)
	    out[k] |= in[i + j] << (2 * j);
      k++;
   }
}

/*
 * n: length of "in" buffer
 *
 * out must be a pointer of length sizeof(char) * n * PACK_DENSITY
 */
void decode(unsigned char *out, const unsigned char *in, const int n)
{
   int i, k;
   unsigned char tmp;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      out[k++] = (tmp & MASK0); 
      out[k++] = (tmp & MASK1) >> 2; 
      out[k++] = (tmp & MASK2) >> 4; 
      out[k]   = (tmp & MASK3) >> 6; 
   }
}

/* 
 *                   plink               sparsnp
 * homozyous:        00 => numeric 0     00 => numeric 0
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * homozygous other: 11 => numeric 3     10 => numeric 2
 * missing:          01 => numeric 1     11 => numeric 3
 *
 * Note that this is reverse order to what the table in
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says, since 
 * the bytes in plink are read backwards HGFEDCBA, not GHEFCDAB
 *
 */
void decode_plink(unsigned char *out, const unsigned char *in, const int n)
{
   int i, k;
   unsigned char tmp;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      
      out[k] = (tmp & MASK0); 
      out[k] = (out[k] == 1) ? 3 : ((out[k] == 3) ? 2 : (out[k] == 2 ? 1 : 0));
      k++;

      out[k] = (tmp & MASK1) >> 2; 
      out[k] = (out[k] == 1) ? 3 : ((out[k] == 3) ? 2 : (out[k] == 2 ? 1 : 0));
      k++;

      out[k] = (tmp & MASK2) >> 4; 
      out[k] = (out[k] == 1) ? 3 : ((out[k] == 3) ? 2 : (out[k] == 2 ? 1 : 0));
      k++;

      out[k] = (tmp & MASK3) >> 6; 
      out[k] = (out[k] == 1) ? 3 : ((out[k] == 3) ? 2 : (out[k] == 2 ? 1 : 0));
   }
}

