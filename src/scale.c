#include "cd.h"

int scale(gmatrix *g, char* filename)
{
   FILE *fileout = NULL;
   
   FOPENTEST(fileout, filename, "w")
      
   return SUCCESS;
}

int main(int argc, char* argv[])
{
   int i, n, p;
   char *filename_in = NULL;
   char *filename_out = NULL;
   gmatrix g;
   short rowmajor = FALSE;
   

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

   }

   if(!gmatrix_init(&g, FALSE, FALSE, rowmajor, filename_in, NULL, NULL, n, p))
      return EXIT_FAILURE;

   scale(&g, filename_out);

   gmatrix_free(&g);

   return EXIT_SUCCESS;
}

