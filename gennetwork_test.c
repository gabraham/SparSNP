/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <stdlib.h>
#include <stdio.h>

#include "gennetwork.h"
#include "common.h"
#include "util.h"

int main(void)
{
   double Ctrue[50] = {
     0.0254602346530758, 0.242074066564266, 0, 0.0918537764912165,
     0, 0, 0.121118658982894, 0, 0, 0, 0.0254602346530758, 0, 0.363657463583651,
     0, 0.016685044285222, 0, 0, 0.0100511989329517, 0, 0, 0, 0.242074066564266,
     -0.363657463583651, 0, 0, 0.0444763436125852, 0, 0, 0.107528132187069,
     0, 0, 0, 0, -0.0918537764912165, 0.016685044285222, -0.0444763436125852,
     0, 0, 0, 3.04866694516442e-05, 0, 0, 0, 0, 0, 0, 0.121118658982894,
     -0.0100511989329517, 0.107528132187069, -3.04866694516442e-05
   };

   /* one-based indexing from R */
   int pairstrue[20] = {
      1L, 1L, 2L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 2L, 3L, 3L,
      4L, 4L, 4L, 5L, 5L, 5L, 5L
   };

   /* one-based indexing from R */
   int edgestrue[20] = {
      1L, 2L, 4L, 7L, 1L, 3L, 5L, 8L, 2L, 3L, 6L, 9L, 4L,
      5L, 6L, 10L, 7L, 8L, 9L, 10L
   };

   /* column-major matrix */
   double Y[100] = {
   1.42, 0.55, 0.84, -0.33, 1.37, -1.78, 1.56, -1.75, 
   -1.49, 1.63, 0.47, 2.94, -0.31, -0.5, -1.37, 1.72, -0.2, -1.09, 
   2.08, -1.36, -0.23, -0.36, -1.92, -0.72, -0.67, 3.25, -1.68, 
   -0.06, -0.16, 0.59, 0.5, 1.58, -0.18, -0.6, -1.44, 0.55, -1.01, 
   1.19, 0.37, -0.66, 0.24, 0.35, -1.93, -0.89, -1.45, 1.55, 0.25, 
   0.63, -1.07, -0.93, 0.23, -0.91, -0.51, 1, 0.43, 0.29, -0.36, 
   0.06, 1.56, 0.37, -0.82, 0.35, 1.42, 1.05, -1.04, 0.62, -0.07, 
   -0.16, -0.22, -1.26, -0.58, -0.05, 0.27, -0.62, -1.09, 1.13, 
   -0.88, -1.51, 0.84, 0.64, -0.34, 0.79, -0.59, -0.03, -0.84, 0.15, 
   1.63, -0.11, 1.06, 0.7, -1.51, 1.02, 0.74, -0.37, -0.31, -0.74, 
   -0.16, 0.47, -0.83, 1.35};
   int N = 20, K = 5;
   int nE = K * (K - 1) / 2;
   int i;
   double s;

   double *C = NULL;
   int *pairs, *edges;

   CALLOCTEST(C, nE * K, sizeof(double));
   CALLOCTEST(pairs, nE * 2, sizeof(int));
   CALLOCTEST(edges, (K - 1) * K, sizeof(int));

   gennetwork(Y, N, K, 0.0, CORTYPE_ABS, C, pairs, edges);

   writematrixf(Y, N, K, "Y.txt");
   writematrixf(C, nE, K, "C.txt");
   writematrixl(pairs, nE, 2, "pairs.txt");
   writematrixl(edges, K - 1, K, "edges.txt");

   s = 0;
   for(i = 0 ; i < nE * K ; i++)
      s += (Ctrue[i] - C[i]) * (Ctrue[i] - C[i]);
   s /= nE * K;
   printf("comparing C and known C, mean square error: %.6f\n", s);

   /* correct for one-based indexing */
   s = 0;
   for(i = 0 ; i < nE * 2 ; i++)
      s += (pairstrue[i] - 1 - pairs[i]) * (pairstrue[i] - 1 - pairs[i]);
   s /= nE * 2;
   printf("comparing pairs and known pairs, mean square error: %.6f\n", s);

   /* correct for one-based indexing */
   s = 0;
   for(i = 0 ; i < K * (K - 1) ; i++)
      s += (edgestrue[i] - 1 - edges[i]) * (edgestrue[i] - 1 - edges[i]);
   s /= nE * 2;
   printf("comparing edges and known edges, mean square error: %.6f\n", s);



   FREENULL(C);
   FREENULL(pairs);
   FREENULL(edges);

   return EXIT_SUCCESS;
}

