
#define 
typedef struct sparsematrix {
   int rows;
   int cols;
   int size;
   int occupied;
   double threshold;
   struct bucket *buckets = NULL;
} sparsematrix;

/*typedef struct sparsevector {
   int size;
   int occupied;
   double threshold;
   struct bucket *buckets = NULL;
} sparsevector;*/

typedef struct bucket {
   short empty;
   int key;
   double value;
} bucket;

/*int sparsevector_set(sparsevector *sv, int i, double x)
{
}

double sparsevector_get(sparsevector *sv, int i)
{
}*/

int bucket_set(bucket *b, int key, double value)
{
   b->key = key;
   b->value = value;
   b->empty = FALSE;
}

int sparsematrix_init(sparsematrix *m, int size, double threshold)
{
   int i;
   m->size = size;
   m->threshold = threshold;
   MALLOCTEST(m->buckets, size * sizeof(bucket));
   for(i = 0 ; i < size ; i++)
      m->buckets[i]->empty = TRUE;

   return SUCCESS;
}

int sparsematrix_free(sparsematrix *m)
{
   if(m->buckets)
   {
      free(m->buckets);
      m->buckets = NULL;
   }
}

int sparsematrix_set(sparsematrix *m, int i, int j, double x)
{
   int key = i * m->cols + j;
   int h = hash(key);

   /* linear probing */
   while(!m->buckets[h]->empty && m->key != key)
   {
      h++;
      if(h == m->size)
	 sparsematrix_grow(m);
   }

   bucket_set(m->buckets + h, key, x);

   m->occupied++;
   return SUCCESS;
}

short sparsematrix_has(sparsematrix *m, int i, int j);
{
   int key = i * m->cols + j;
   int h = hash(key);

   if(m->buckets[h]->empty)
      return FALSE;

   while(!m->buckets[h]->empty)
   {
      h++;
   }

   return FALSE;
}

double sparsematrix_get(sparsematrix *m, int i, int j)
{
   int key = i * m->cols + j;
   int h = hash(key);

   if(m->buckets[h]->empty)
      return NaN;

   /* linear probing */
   while(m->buckets[h]->empty || m->key != key)
      h++;

   return m->buckets + h, key, x);
}

int sparsematrix_grow(sparsematrix *m)
{
   m->size *= 2;
   REALLOCTEST(m->buckets, m->buckets, m->size * sizeof(bucket))
   return SUCCESS;
}

int sparsematrix_free(sparsematrix *m)
{
}

