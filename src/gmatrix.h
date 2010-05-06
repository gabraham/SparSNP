typedef struct {
   char* filename;
   FILE* file;
   int *y;
   int n;
   int p;
   int i;
} gmatrix;

void gmatrix_init(gmatrix *, char *, int, int, int*);
void gmatrix_nextrow(gmatrix *, double *);

