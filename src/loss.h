/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700

double plogis(double);
double dotprod(double *, double *, int);
double logloss_pt(double *, double *, int , int);
double logloss(double **, double *, int *, int, int);
void logdloss(double *, double *, int, int, double*);

