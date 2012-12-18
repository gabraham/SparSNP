/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#define THIN_WINDOW_SIZE 50
#define THIN_STEP_SIZE 5
#define THIN_COR_MAX 0.90 /* |r|, not r^2 */

void copyshrink(double *x, double *y, int n, int p, int *active, int m);
int thin(double *x, int n, int p, int *active,
      double cor_max, int windowsize, int stepsize);


