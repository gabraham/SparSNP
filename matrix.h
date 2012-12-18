/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

void crossprod(double *x, double *y, double *z, int m, int n, int p);
void invert2x2(double *y, double *x);
void wcrossprod(const double *x, const double *y, const double *w, double *z,
      const int m, const int n, const int p);
void sqmvprod(const double *x, const double *y, double *z, int m);
void copyshrink(double *x, double *y, int n, int p, int *active, int m);
void copyshrinkrange(double *x, double *y, int n, int p, int from, int to);
void printmatrix(double *x, int n, int p);
void printmatrix0(double *x, int n, int p);
int cov(double *x, double *S, int n, int p);
int cov2(double *x, double *S, int n, int p);
void cov2cor(double *S, double *P, int p);
int covmiss(double *x, double *S, int n, int p, int *good);

