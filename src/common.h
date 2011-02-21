#include <stdlib.h>
 
#define VERBOSE 1

#define sign(x) ((x > 0) - (x < 0))
#define TRUE 1
#define FALSE 0

#define SUCCESS 1
#define FAILURE 0

#define MAX_STR_LEN 101

/* Value below which std dev is considered zero */
#define SDTHRESH 1e-10

/* Value below which is considered zero, for use in
 * absolute convergence test */
#define ZERO_THRESH 1e-7

#define ONE 1.0

#define CDFAILURE -1

/* The size for each datum in the binary input file */
#ifndef DTYPE
#define DTYPE unsigned char
#endif

#ifndef DTYPE_DEF
typedef DTYPE dtype;
#define DTYPE_DEF 1
#endif

#define MODEL_LINEAR 1
#define MODEL_PCOR 2
#define MODEL_LOGISTIC 3
#define MODEL_SQRHINGE 4

#define MODEL_NAME_LINEAR "linear"
#define MODEL_NAME_PCOR "pcor"
#define MODEL_NAME_LOGISTIC "logistic"
#define MODEL_NAME_SQRHINGE "sqrhinge"

#define MODE_TRAIN 0
#define MODE_PREDICT 1

/* format of binary input */
#define BINFORMAT_BIN 1
#define BINFORMAT_PLINK 2
	  
/* number of byts in plink BED header */
#define PLINK_HEADER_SIZE 3

#define FMAX(a, b) (a < b ? b : a) 
#define FMIN(a, b) (a < b ? a : b)

/* Macros with built in error checking */

size_t retval;

#define MALLOCTEST(x, size) \
if(!(x = malloc(size))) { \
fprintf(stderr, "can't malloc\n"); \
return FAILURE; \
}

#define MALLOCTEST2(x, size) \
if(!(x = malloc(size))) { \
fprintf(stderr, "can't malloc\n"); \
return EXIT_FAILURE; \
}

#define CALLOCTEST(x, count, size) \
if(!(x = calloc(count, size))) { \
fprintf(stderr, "can't calloc\n"); \
return FAILURE; \
}

#define CALLOCTEST2(x, count, size) \
if(!(x = calloc(count, size))) { \
fprintf(stderr, "can't calloc\n"); \
return EXIT_FAILURE; \
}

#define REALLOCTEST(x, y, size) \
if(!(x = realloc(y, size))) { \
fprintf(stderr, "can't realloc\n"); \
return FAILURE; \
}

#define REALLOCTEST2(x, y, size) \
if(!(x = realloc(y, size))) { \
fprintf(stderr, "can't realloc\n"); \
return EXIT_FAILURE; \
}

#define FOPENTEST(x, filename, mode) \
if(!(x = fopen(filename, mode))) { \
fprintf(stderr, "can't open file %s\n", filename); \
return FAILURE; \
}

#define FSEEKOTEST(x, offset, whence) \
if(fseeko(x, offset, whence) != 0) { \
fprintf(stderr, "can't seek offset %lld\n", \
(unsigned long long)offset); \
return FAILURE; \
}

#define FREADTEST(x, size, count, stream) \
if((retval = fread(x, size, count, stream)) < count) { \
fprintf(stderr, "read fewer items (%lu) than expected (%lu)\n", \
retval, (unsigned long)count); \
fflush(stderr); return FAILURE; \
}

#define FWRITETEST(x, size, count, stream) \
if((retval = fwrite(x, size, count, stream)) < count) { \
fprintf(stderr, "wrote fewer items (%lu) than expected (%lu)\n", \
retval, (unsigned long)count); \
return FAILURE; \
}

#define FREENULL(x) \
if(x) { \
free(x); \
x = NULL; \
}

int strcmp2(const char*, const char*);

double soft_threshold(double beta, double gamma);

