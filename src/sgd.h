void writebeta(char*, double*, int);
double sgd_gmatrix(gmatrix *, double,
      int, double *, double, double,
      double, int, int *, double);
void predict_logloss(gmatrix *, double *, double *, int *);
double predict_logloss_pt(sample *, double *, double *, double *, int);

