
#define sign(x) ((x > 0) - (x < 0))
#define TRUE 1
#define FALSE 0

#define SUCCESS 1
#define FAILURE 0

/* How to treat the x input, as discrete or continuous */
#ifdef DISCRETE
#define type "discrete"
#define dtype char
#define ONE 1
#else
#define dtype double
#define type "continuous"
#define ONE 1.0
#endif

/* The size for each datum in the binary input file */
#define intype char


int strcmp2(const char*, const char*);


