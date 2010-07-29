
#define VERBOSE 1

#define sign(x) ((x > 0) - (x < 0))
#define TRUE 1
#define FALSE 0

#define SUCCESS 1
#define FAILURE 0

#define MAX_STR_LEN 100

/* Value below which std dev is considered zero */
#define SDTHRESH 1e-10

/* Value below which is considered zero, for use in
 * absolute convergence test */
#define ZERO_THRESH 1e-12

#define ONE 1.0

/* The size for each datum in the binary input file */
#ifndef DTYPE
#define DTYPE char
#endif

#ifndef DTYPE_DEF
typedef DTYPE dtype;
#define DTYPE_DEF 1
#endif


/* Macros with built in error checking */

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

#define FOPENTEST(x, filename, mode) \
if(!(x = fopen(filename, mode))) { \
fprintf(stderr, "can't open file %s\n", filename); \
return FAILURE; \
}

#define FSEEKTEST(x, offset, whence) \
if(fseek(x, offset, whence) != 0) { \
fprintf(stderr, "can't seek offset %ld\n", offset); \
return FAILURE; \
}

#define FREADTEST(x, size, count, stream) \
if(fread(x, size, count, stream) < count) { \
fprintf(stderr, "read fewer items than expected (%d)\n", count); \
return FAILURE; \
}

#define FWRITETEST(x, size, count, stream) \
if(fwrite(x, size, count, stream) < count) { \
fprintf(stderr, "wrote fewer items than expected (%d)\n", count); \
return FAILURE; \
}

int strcmp2(const char*, const char*);

double soft_threshold(double beta, double gamma);

