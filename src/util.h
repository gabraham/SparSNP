int writevectorf(char* file, double* beta, int p);
int writevectorl(char* file, int* beta, int p);
int writematrixf(double **x, int n, int p, char* file);
int load_beta(double *beta, char *filename, int p);
void scale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p);
void unscale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p);
int writebinvectorf(char* file, double* x, int p);
int writebinvectorl(char* file, int* x, int p);
int readvectorl(char *file, int *x, int n);
