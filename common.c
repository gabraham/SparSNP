/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
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

