typedef double (*loss_pt)(dtype *, double *, dtype, int);
typedef double (*predict_pt)(sample *, double *, double *, double *, int);
typedef void (*dloss)(dtype *, double *, dtype, int, double *);

void writebeta(char*, double*, int);
void predict_logloss(gmatrix *, double *, double *, int *);
double predict_logloss_pt(sample *, double *, double *, double *, int);
double sgd_gmatrix(gmatrix *g,
   dloss, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);


