/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

void tcrossprod(double *A, double *B, int *m, int *k, int *n, double *C);
int pseudoinverse(double *Aorig, int *m, int *n, double *P);
   

