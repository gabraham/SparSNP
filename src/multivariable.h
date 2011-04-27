#define NEWTON_THRESH 1e-3
#define NEWTON_THRESH_MAX 30
#define NEWTON_DEVIANCE_MIN 1e-2

#define NEWTON_SUCCESS 1

/* didn't converge within predefined iterations */
#define NEWTON_ERR_NO_CONVERGENCE 2 

/* converged but to a very large value */
#define NEWTON_ERR_DIVERGENCE 3 

int newton(double *x, double *y, double *beta, double *invhessian,
      int n, int p, double lambda2, int verbose);

int make_hessian(double *hessian, double *x,
      double beta_intercept, double beta, int n);

int multivariable_newton(Opt *opt, gmatrix *g, int nums1,
      int *pselected, int *numselected, int *rets);

int multivariable_lasso(Opt *opt, gmatrix *g, int threshind);
int make_lambda1path(Opt *opt, gmatrix *g, int threshind);

