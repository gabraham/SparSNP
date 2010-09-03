
#define PACK_DENSITY 4

void encode(unsigned char *out, const unsigned char *in,
      const int n);
void decode(unsigned char *out, const unsigned char *in,
      const int n);
void decode_plink(unsigned char *out, const unsigned char *in,
      const int n);

