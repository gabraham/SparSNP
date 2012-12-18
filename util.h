/* 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */
 
int writevectorf(char* file, double* beta, int p);
int writevectorl(char* file, int* beta, int p);
int writematrixf(double *x, int n, int p, char* file);
int writematrixl(int *x, int n, int p, char* file);

int load_beta(double *beta, char *filename, int p);
void scale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p);
void unscale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p, int K);
int writebinvectorf(char* file, double* x, int p);
int writebinvectorl(char* file, int* x, int p);
int readvectorl(char *file, int *x, int n);
int write_beta_sparse(char* file, double* beta, int p, int K);
int load_beta_sparse(double *beta, char *filename, int p);
