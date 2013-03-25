/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <math.h>
#include <stdio.h>

#include "gennetwork.h"
#include "common.h"
#include "matrix.h"
#include "util.h"

/*
 * Creates the nE by nV (num edges by num vertices/tasks) fusion penalty
 * matrix
 * 
 * pairs is an nE * 2 matrix of which tasks are in which pair
 */
int gennetwork(double *y, int n, int K,
   double corthresh, int cortype, double *C, int *pairs, int *edges)
{
   int i, j, k, e,
      nE = K * (K - 1) / 2;
   double *S = NULL,
	  *R = NULL,
	  *W = NULL;
   int *eFrom = NULL,
       *eTo = NULL;

   CALLOCTEST(S, K * K, sizeof(double));
   CALLOCTEST(R, K * K, sizeof(double));
   CALLOCTEST(W, K * K, sizeof(double));
   CALLOCTEST(eFrom, nE, sizeof(int));
   CALLOCTEST(eTo, nE, sizeof(int));

   cov2(y, S, n, K);
   cov2cor(S, R, K);

#ifdef DEBUG
   if(!writematrixf(S, K, K, "S.txt"))
      return FAILURE;
   if(!writematrixf(R, K, K, "R.txt"))
      return FAILURE;
#endif

   for(i = 0 ; i < K ; i++)
   {
      for(j = 0 ; j < K ; j++)
      {
	 k = i * K + j;
         R[k] = fabs(R[k]) < corthresh ? 0 : R[k];
   
         if(cortype == CORTYPE_IND)
	    W[k] = fabs(R[k]) > corthresh;
	 else if(cortype == CORTYPE_ABS)
	    W[k] = fabs(R[k]);
	 else if(cortype == CORTYPE_SQR)
	    W[k] = R[k] * R[k];
      }
   }

   e = 0;
   for(j = 0 ; j < K ; j++)
   {
      for(i = 0 ; i < K ; i++)
      {
	 if(i < j && R[i * K + j] != 0)
	 {
	    eFrom[e] = i;
	    eTo[e] = j;
	    pairs[e] = i;
	    pairs[e + nE] = j;
	    e++;
	 }
      }
   }

   for(e = 0 ; e < nE ; e++)
   {
      i = eFrom[e];
      j = eTo[e];

      C[e + nE * i] = W[i * K + j];
      C[e + nE * j] = -W[i * K + j] * sign(R[i * K + j]);
   }

   /* make the (K-1) by K edges matrix */
   for(j = 0 ; j < K ; j++)
   {
      e = 0;
      for(i = 0 ; i < nE ; i++)
      {
	 if(C[i + nE * j] != 0)
	 {
	    edges[e + (K - 1) * j] = i;
	    e++;
	 }
      }
   }

   FREENULL(S);
   FREENULL(R);
   FREENULL(W);
   FREENULL(eFrom);
   FREENULL(eTo);

   return SUCCESS;
}

