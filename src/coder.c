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
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 * 
 */
void decode_plink(unsigned char *out, const unsigned char *in, const int n)
{
   int i, k;
   unsigned char tmp, geno;
   int a1, a2;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      
      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */
      geno = (tmp & MASK0);
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK1) >> 2; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK2) >> 4; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK3) >> 6; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}

