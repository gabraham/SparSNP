/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>    
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/param.h>

int main(int argc, char *argv[])
{
   char *symlinkpath = NULL;
   char actualpath [PATH_MAX];
   char *ptr;

   if(argc > 1)
   {
      symlinkpath = argv[1];
      if((ptr = realpath(symlinkpath, actualpath)))
      {
	 printf("%s\n", ptr);
	 return EXIT_SUCCESS;
      }
      fprintf(stderr, "realpath: %s when reading '%s'\n",
	 strerror(errno), symlinkpath);
   }
   else
      fprintf(stderr, "realpath: missing file name\n");
   
   return EXIT_FAILURE;
}

