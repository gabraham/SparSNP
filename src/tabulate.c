#include "cd.h"
#include "util.h"

int main(int argc, char* argv[])
{
   int i, n = 0, p = 0;
   int xlevels = 3, ylevels = 2;
   char *filename_in = NULL;
   char *filename_out = NULL;
   short inmemory = FALSE;
   short doscale = FALSE;
   short verbose = FALSE;
   gmatrix g;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-fin"))
      {
	 i++;
	 filename_in = argv[i];
      }
      else if(strcmp2(argv[i], "-fout"))
      {
	 i++;
	 filename_out = argv[i];
      }
      else if(strcmp2(argv[i], "-n"))
      {
	 i++;
	 n = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-p"))
      {
	 i++;
	 p = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-inmemory"))
      {
	 inmemory = TRUE;
      }
      else if(strcmp2(argv[i], "-scale"))
      {
	 doscale = TRUE;
      }
      else if(strcmp2(argv[i], "-xlev"))
      {
	 i++;
	 xlevels = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-ylev"))
      {
	 i++;
	 ylevels = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-v"))
      {
	 verbose = TRUE;
      }
   }

   if(filename_in == NULL || filename_out == NULL || n == 0
      || p == 0 || xlevels == 0 || ylevels == 0)
   {
      printf("tabulate: -fin <filein> -fout <fileout> -n #n -p #p\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_in, n, p, inmemory, FALSE, NULL))
      return EXIT_FAILURE;

   tabulate(&g, filename_out, doscale, xlevels, ylevels, verbose);

   /*writevectorf("mean.csv", g.mean, p + 1);
   writevectorf("sd.csv", g.sd, p + 1);*/

   gmatrix_free(&g);

   return EXIT_SUCCESS;
}

