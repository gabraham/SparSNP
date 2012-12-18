/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#define CORTYPE_IND 0
#define CORTYPE_ABS 1
#define CORTYPE_SQR 2

int gennetwork(double *y, int n, int K,
   double corthresh, int cortype, double *C);


