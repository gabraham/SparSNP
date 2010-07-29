#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "gmatrix.h"

#define MAXLP 20

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
typedef double (*phi1)(double);
typedef double (*phi2)(double);
typedef double (*inv)(double);
typedef double (*step)(sample *s, double *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func);

typedef struct Opt {
   int maxepochs;
   double lambda1;
   double lambda2;
   double threshold;
   double l1minratio;
   int nlambda1;
   double trunc;
   char *model;
   char *betafile;
   loss_pt loss_pt_func;
   phi1 phi1_func;
   phi2 phi2_func;
   inv inv_func;
   short nofit;
   int n, p;
   short warmrestarts;
   char *filename;
   double *lambda1path;
   double lambda1max;
   double lambda1min;
   int verbose;
   int cv;
   long seed;
   int nzmax;
   int ntrain;
   int *trainf;
   char *subsetfile;
   char *lambda1pathfile;
   step step_func;
   short inmemory;
   short tabulate;
   char *scalefile;
} Opt;

int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
      inv inv_func,
      step step_func,
      int maxepoch, double *beta, double *lp,
      double lambda1, double lambda2,
      double threshold, int verbose, int *trainf, double trunc);

/*int cd_gmatrix2(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,
      int maxepoch,
      double *beta,
      double lambda1);*/

double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      inv inv_func,
      step step_func
);

int cvsplit(Opt *opt);
void opt_free(Opt *opt);
void opt_defaults(Opt *opt);
int opt_parse(int argc, char* argv[], Opt* opt);
int make_lambda1path(Opt *opt, gmatrix *g);
int run(Opt *opt, gmatrix *g);

double step_regular(sample *s, double *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func);

double step_grouped(sample *s, double *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func);

