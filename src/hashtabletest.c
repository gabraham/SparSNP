#include <stdlib.h>
#include <stdio.h>
#include "hashtable.h"

int main()
{
   int i;
   double *values, *keys;
   double val;
   hashtable ht;
   int N = 12000;

   if(!hashtable_init(&ht))
      return EXIT_FAILURE;

   MALLOCTEST2(values, sizeof(double) * N)
   MALLOCTEST2(keys, sizeof(double) * N)

   for(i = 0 ; i < N ; i++)
   {
      keys[i] = (double)i;
      values[i] = 2.0 * i;
   }

   for(i = 0 ; i < N ; i++)
      if(!hashtable_put(&ht, keys[i], values[i]))
	 return EXIT_FAILURE;

   for(i = 0 ; i < N ; i++)
   {
      val = hashtable_get(&ht, keys[i]);
      if(val != values[i])
      {
	 printf("mismatch %.3f: is %.3f, should be %.3f\n",
	       keys[i], val, values[i]);
      }
      else
	 printf("match %.3f: %.3f\n", keys[i], val);
   }

   free(values);
   free(keys);

   hashtable_free(&ht);

   return 0;
}

