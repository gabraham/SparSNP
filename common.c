/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <string.h>
#include <math.h>
#include "common.h"

int strcmp2(const char* s1, const char* s2)
{
   return strcmp(s1, s2) == 0;
}

double soft_threshold(double beta, double gamma)
{
   return sign(beta) * fmax(fabs(beta) - gamma, 0);
}

