#include "common.h"

unsigned int MurmurHash2(const void *key, int len, unsigned int seed);

typedef struct bucket {
   double key;
   double value;
   unsigned int active;
   struct bucket *next;
} bucket;

typedef struct hashtable {
   unsigned int size;
   unsigned int active;
   bucket *buckets;
} hashtable;

int hashtable_init(hashtable *ht);
void hashtable_free(hashtable *ht);
int hashtable_put(hashtable *ht, double key, double value);
double hashtable_get(hashtable *ht, double key);

