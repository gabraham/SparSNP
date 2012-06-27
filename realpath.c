
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>    

int
main(int argc, char *argv[])
{
    char *symlinkpath = argv[1];
    char actualpath [PATH_MAX];
    char *ptr;
    ptr = realpath(symlinkpath, actualpath);
    printf("%s\n", ptr);
    return EXIT_SUCCESS;
}

