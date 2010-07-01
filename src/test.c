#include <stdlib.h>
#include <stdio.h>

int main()
{
   int i = 0, j;
   char *x = malloc(sizeof(char) * 20);
   FILE* f = fopen("x.bin.t", "rb");
   
   while(fread(x, sizeof(char), 20, f))
   {
      printf("%d: ", i);
      for(j = 0 ; j < 20 ; j++)
	 printf("%.1f ", (double)x[j]);
      printf("\n");
      i++;
   }

   printf("---\n");

   return EXIT_SUCCESS;
}

