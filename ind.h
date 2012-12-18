/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

int ind_read(char *file, int *folds, const int n, const int nfolds);
int ind_write(char *file, int *folds, int n, int nfolds);
int ind_getfolds(char *file);

