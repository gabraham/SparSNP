
#define PACK_DENSITY 4

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

typedef struct mapping {
   int **map;
} mapping;

void encode(unsigned char *out, const unsigned char *in,
      const int n);
void decode(unsigned char *out, const unsigned char *in,
      const int n);
void decode_plink(unsigned char *out, const unsigned char *in,
      const int n);
int mapping_init(mapping *m);
void mapping_free(mapping *m);

void decode_plink_mapping(mapping *m, unsigned char *out,
      const unsigned char *in, const int n);

