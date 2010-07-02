
#define sign(x) ((x > 0) - (x < 0))
#define TRUE 1
#define FALSE 0

#define SUCCESS 1
#define FAILURE 0

#define SDTHRESH 1e-10

/* How to treat the x input, as discrete or continuous */
/*#ifdef DISCRETE
#define type "discrete"
#define dtype char
#define ONE 1
#else
#define dtype float
#define type "continuous"
#define ONE 1.0
#endif */

#define dtype double
#define type "continuous"
#define ONE 1.0

/* The size for each datum in the binary input file */
#ifndef intype
#define intype char
#endif

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

#define FREADTEST(x, size, count, stream) \
if(fread(x, size, count, stream) < count) { \
fprintf(stderr, "read fewer bytes than expected (%d)\n", count); \
return FAILURE; \
}

int strcmp2(const char*, const char*);

double soft_threshold(double beta, double gamma);

