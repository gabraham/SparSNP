#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "gmatrix.h"

#define MAXLP 7

#define MAX_SHOW_NOTCONV 100

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
typedef double (*phi1)(double);
typedef double (*phi2)(double);
typedef double (*inv)(double);
typedef double (*step)(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func);

typedef double (*predict)(double x);

typedef struct Opt {
   short mode;
   short model;
   int maxepochs;
   int maxiters;
   double lambda1;
   double lambda2;
   double threshold;
   double l1minratio;
   int nlambda1;
   double trunc;
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
   int nfolds;
   long seed;
   int nzmax;
   int ntrain;
   int *trainf;
   char *subsetfile;
   char *lambda1pathfile;
   step step_func;
   short inmemory;
   char *scalefile;
   short yformat;
   predict predict_func;
   char **beta_files;
   int n_beta_files;
   char *predict_file;
   short encoded;
   short binformat;
   char *folds_ind_file;
} Opt;

int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      step step_func,
      const int maxepochs, const int maxiters,
      const double lambda1, const double lambda2,
      const double threshold, const int verbose,
      const double trunc);

double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      inv inv_func,
      step step_func
);

int cvsplit(Opt *opt);
void opt_free(Opt *opt);
int opt_defaults(Opt *opt);
int opt_parse(int argc, char* argv[], Opt* opt);
int make_lambda1path(Opt *opt, gmatrix *g);
int run(Opt *opt, gmatrix *g);
void zero_model(gmatrix *g);

double step_generic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func);

double step_regular_linear(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func);

double step_regular_logistic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func);

double step_regular_sqrhinge(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func);



