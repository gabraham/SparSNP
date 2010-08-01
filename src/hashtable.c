#include "cd.h"
#include "util.h"
#include "hashtable.h"

#ifndef HASH
#define HASH MurmurHash2
#define SEED 984634157 /* arbitrary */
#endif

#define HASH_SIZE 1048576
#define MOD_SHIFT 12 /* 2**32 / 2**20 = 2**12 bits*/
#define EXP_ERROR 1e-6

int hashtable_init(hashtable *ht)
{
   unsigned int i;
   ht->size = HASH_SIZE;
   ht->active = 0;

   MALLOCTEST(ht->buckets, sizeof(bucket) * ht->size)

   for(i = 0 ; i < ht->size ; i++)
   {
      ht->buckets[i].active = FALSE;
      ht->buckets[i].next = NULL;
      ht->buckets[i].key = NAN;
   }

   return SUCCESS;
}

void hashtable_free(hashtable *ht)
{
   unsigned int i;
   bucket *bk1, *bk2;

   if(ht->buckets)
   {
      /* free the chain except for the first bucket */
      for(i = 0 ; i < ht->size ; i++)
      {
	 bk1 = ht->buckets[i].next;
	 while(bk1)
	 {  
	    bk2 = bk1->next;
	    free(bk1);
	    bk1 = bk2;
	 }
      }

      /* free the first buckets */
      free(ht->buckets);
   }
   ht->buckets = NULL;
   ht->active = 0;
}

/* If a key already exists, we don't overwrite it */
int hashtable_put(hashtable *ht, double key, double value)
{
   bucket *bk = NULL, *bk_prev = NULL;
   unsigned int hval = HASH(&key, sizeof(double), SEED);

   /* traverse chain */
   bk = &ht->buckets[hval];
   while(bk && bk->active && fabs(bk->key - key) > EXP_ERROR)
   {
      bk_prev = bk;
      bk = bk->next;
   }

   if(!bk)
   {
      MALLOCTEST(bk_prev->next, sizeof(bucket));
      bk = bk_prev->next;
      bk->next = NULL;
   }

   bk->key = key;
   bk->value = value;
   bk->active = TRUE;
   ht->active++;

   return SUCCESS;
}

double hashtable_get(hashtable *ht, double key)
{
   bucket *bk;
   unsigned int hval = HASH(&key, sizeof(double), SEED);

   /* traverse chain */
   bk = &ht->buckets[hval];
   while(bk && bk->active && fabs(bk->key - key) > EXP_ERROR)
   {
      bk = bk->next;
   }

   if(bk && bk->active)
      return bk->value;

   return NAN;
}

//-----------------------------------------------------------------------------
// MurmurHash2, by Austin Appleby

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed )
{
   // 'm' and 'r' are mixing constants generated offline.
   // They're not really 'magic', they just happen to work well.

   const unsigned int m = 0x5bd1e995;
   const int r = 24;

   // Initialize the hash to a 'random' value

   unsigned int h = seed ^ len;

   // Mix 4 bytes at a time into the hash

   const unsigned char * data = (const unsigned char *)key;

   while(len >= 4)
   {
      unsigned int k = *(unsigned int *)data;

      k *= m; 
      k ^= k >> r; 
      k *= m; 

      h *= m; 
      h ^= k;

      data += 4;
      len -= 4;
   }

   // Handle the last few bytes of the input array

   switch(len)
   {
      case 3: h ^= data[2] << 16;
      case 2: h ^= data[1] << 8;
      case 1: h ^= data[0];
	      h *= m;
   };

   // Do a few final mixes of the hash to ensure the last few
   // bytes are well-incorporated.

   h ^= h >> 13;
   h *= m;
   h ^= h >> 15;

   /*return h;*/
   return h >> MOD_SHIFT;
} 


