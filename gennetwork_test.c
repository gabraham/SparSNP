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

   double *C = NULL;

   CALLOCTEST(C, 5 * 4 / 2 * 5, sizeof(double));

   gennetwork(Y, 20, 5, 0.0, CORTYPE_ABS, C);

   writematrixf(C, 10, 5, "C.txt");


   FREENULL(C);

   return EXIT_SUCCESS;
}

