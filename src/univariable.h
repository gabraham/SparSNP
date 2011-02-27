#define OPTIONS_CALLER univariable

#define NR_THRESH 1e-3
#define NR_THRESH_MAX 30
#define NR_DEVIANCE_MIN 1e-2

/* didn't converge within predefined iterations */
#define NR_ERR_NO_CONVERGENCE 2 

/* converged but to a very large value */
#define NR_ERR_DIVERGENCE 3 

int run_train_nr(Opt *opt, gmatrix *g, int nums1,
      int *pselected, int *numselected, int *rets);

int run_train_lasso(Opt *opt, gmatrix *g, int nums1,
      int *pselected, int *numselected, int *rets);

