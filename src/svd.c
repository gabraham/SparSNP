#include "svd.h"

/*
 * This declaration is necessary since /usr/lib/clapack.h on Linux doesn't
 * have it
 */
int dgesdd(char *jobz, int *m, int *n, double *a, int *lda,
      double *s, double *u, int *ldu, double *vt, int *ldvt,
      double *work, int *lwork, int *iwork)
{
   extern void dgesdd_(const char *JOBZ, const int *M, const int *N,
	 const double *A, const int *LDA, double* S, double* U, const int *LDU,
	 double *VT, int *LDVT, double *WORK, int *LWORK, int *IWORK,
	 int *INFO);

   int info;
   dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt,
	 work, lwork, iwork, &info);
   return info;
}

/*
 * C = A B'
 */
void tcrossprod(double *A, double *B, int *m, int *k, int *n, double *C)
{
   double alpha = 1.0, beta = 0.0;
   int LDA = *m;
   int LDB = *n;
   int LDC = *m;

   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	 *m, *n, *k, alpha, A, LDA, B, LDB, beta, C, LDC);
}

/*
 * Takes the m by n matrix x, drops the columns for which the flags are 0,
 * returns the resulting matrix.
 */
double* dropcolumns(double *x, int m, int n, int* flags)
{
   int i, j, k, num = 0;
   double *x2;

   for(j = 0 ; j < n ; j++)
      num += flags[j];
   
   x2 = Calloc(m * num, double);

   for(i = 0 ; i < m ; i++)
   {
      k = 0;
      for(j = 0 ; j < n ; j++)
      {
	 if(flags[j])
	 {
      	    x2[i + k * m] = x[i + j * m]; 
	    k++;
	 }
      }
   }

   return x2;
}

/*
 *  flags: int vector of length m
 */
double* droprows(double *x, int m, int n, int* flags)
{
   int i, j, k, num = 0;
   double *x2;

   for(i = 0 ; i < m ; i++)
      num += flags[i];

   x2 = Calloc(num * n, double);

   k = 0;
   for(i = 0 ; i < m ; i++)
   {
      if(flags[i])
      {
	 for(j = 0 ; j < n ; j++)
	    x2[k + j * num] = x[i + j * m];
	 k++;
      }
   }

   return x2;
}

/*
 *  SVD: X = U S V'
 *
 *  Pseudoinverse: P = V' S^{-1}' U', i.e., v %*% diag(1/d) %*% t(u)
 *
 *  A: an m * n matrix
 *  P: an n * m matrix
 *
 *  A tolerance is calculated, and eigenvectors with related eigenvalues of
 *  less than the tolerance will be dropped.
 */
int pseudoinverse(double *Aorig, int *m, int *n, double *P)
{
   double *A, *U, *U2, *S, *S2, *VT, *VT2,
	  *WORK, *TMP, optwork = 0, *Sinv, tol, maxs;

   int LDA = *m,
       LDU = *m,
       LDVT,
       LWORK = -1,
       *IWORK,
       INFO,
       i, j, k,
       q,
       *eigens, neigens;
   char JOBZ[] = {'S'};

   A = Calloc((*m) * (*n), double);
   memcpy(A, Aorig, sizeof(double) * (*m) * (*n));

   q = fmin(*m, *n);
   LDVT = q;

   IWORK = Calloc(8 * q, int);
   S = Calloc(q, double);
   eigens = Calloc(q, int);
   U = Calloc(LDU * q, double);
   VT = Calloc(q * (*n), double);

   #ifdef VERBOSE
   if(VERBOSE > 1)
   {
      printf("A\n");
      printmatrix(A, m, n, CblasColMajor);
   }
   #endif

   /* Query to get the optimal work space size */
   INFO = dgesdd(JOBZ, m, n, A, &LDA, S, U, &LDU, VT,
	 &LDVT, &optwork, &LWORK, IWORK);
   if(INFO != 0)
   {
      /*error("Error %d from dgesdd", INFO);*/
      fprintf(stderr, "Error %d from dgesdd", INFO);
      return FAILURE;
   }

   LWORK = (int)optwork;
   WORK = Calloc(LWORK, double);

   /* Now do actual SVD */
   INFO = dgesdd(JOBZ, m, n, A, &LDA, S, U, &LDU, VT,
	 &LDVT, WORK, &LWORK, IWORK);
   Free(IWORK);
   Free(WORK);
   Free(A);
   if(INFO != 0)
   {
      /*error("Error %d from dgesdd", INFO);*/
      fprintf(stderr, "Error %d from dgesdd", INFO);
      return FAILURE;
   }

   /* max(S) */
   maxs = S[0];
   for(i = 1 ; i < q ; i++)
      if(maxs < S[i])
	 maxs = S[i];

   /* Drop eigenvectors with eigenvalues that are close to zero */

   /* Same tolerance heuristic as in corpcor package */
   tol = fmax(*m, *n) * maxs * EPS;

   #ifdef VERBOSE
   if(VERBOSE > 1)
      printf("tol=%5.15f\n", tol);
   #endif

   /* Find which eigenvalues fall above tol (1: above, 0: below) */
   neigens = 0;
   for(i = 0 ; i < q ; i++)
   {
      neigens += eigens[i] = fabs(S[i]) <= tol ? 0 : 1;

      #ifdef VERBOSE
      if(VERBOSE > 1)
	 printf("eigens[%d]=%d\n", i, eigens[i]);
      #endif
   }

   #ifdef VERBOSE
   if(VERBOSE > 1)
      printf("neigens=%d\n", neigens);
   #endif

   /* 
    * Return a matrix of zeros of dimensions n * m
    */
   if(neigens == 0)
   {
      #ifdef VERBOSE
      if(VERBOSE > 1)
	 printf("All eigenvalues are zero, returning zero matrix\n");
      #endif

      for(i = 0 ; i < *m ; i++)
	 for(j = 0 ; j < *n ; j++)
	    P[i + j * (*m)] = 0;

      Free(U);
      Free(VT);
      return SUCCESS;
   }

   U2 = dropcolumns(U, LDU, q, eigens);
   VT2 = droprows(VT, q, *n, eigens);
   Sinv = calloc(neigens * neigens, sizeof(double));
   Free(U);
   Free(VT);

   #ifdef VERBOSE
   if(VERBOSE > 1)
   {
      printf("U\n");
      printmatrix(U, &LDU, &q, CblasColMajor);
      printf("VT\n");
      printmatrix(VT, &q, n, CblasColMajor);

      printf("U2\n");
      printmatrix(U2, &LDU, &neigens, CblasColMajor);
      printf("VT2\n");
      printmatrix(VT2, &neigens, n, CblasColMajor);
   }
   #endif

   /* Sinv = S^{-1}
    *
    * Sinv is a matrix, with the diagonal being the vector S
    * Ignore the tiny eigenvalues
    */
   #ifdef VERBOSE
   if(VERBOSE > 1)
      printf("S\n");
   #endif

   k = 0;
   for(i = 0 ; i < q ; i++)
   {
      #ifdef VERBOSE
      if(VERBOSE > 1)
	 printf("%10.5g ", S[i]);
      #endif

      if(eigens[i])
      {
	 Sinv[k * neigens + k] = 1.0 / S[i];
	 k++;
      }
   }
   Free(S);

   #ifdef VERBOSE
   if(VERBOSE > 1)
   {
      printf("\n");
      printf("Sinv\n");
      printmatrix(Sinv, &neigens, &neigens, CblasColMajor);
   }
   #endif

   /* TMP = VT' Sinv */
   TMP = Calloc((*n) * neigens, double);
   cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	 *n, neigens, neigens, 1.0, VT2, neigens, Sinv, neigens, 0, TMP, *n);
   Free(VT2);
   Free(Sinv);

   #ifdef VERBOSE
   if(VERBOSE > 1)
   {
      printf("TMP\n");
      printmatrix(TMP, n, &neigens, CblasColMajor);
   }
   #endif

   /* P = TMP U' */
   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	 *n, LDU, neigens, 1.0, TMP, *n, U2, LDU, 0, P, *n);
   Free(TMP);
   Free(U2);

   return SUCCESS;
}

