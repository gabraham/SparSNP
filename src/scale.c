#include <math.h>

#include "scale.h"
#include "gmatrix.h"

void scale(gmatrix *g, double *mean, double *sd)
{
   int i, j;
   int n = g->n;
   int p = g->p;
   double delta;
   sample sm;

   sample_init(&sm, p);
 
   /* sd is really the sum of squares, not the SD, but we
    * use the same variable to save memory */

   for(i = 0 ; i < n ; i++)
   {
      gmatrix_nextrow(g, &sm);

     for(j = 0 ; j < p ; j++)
     {
         if(i == 0)
            mean[j] = sd[j] = 0;

         delta = sm.x[j] - mean[j];
         mean[j] += delta / (i + 1);
         sd[j] += delta * (sm.x[j] - mean[j]);
      }
   }

   for(j = 0 ; j < p ; j++)
      sd[j] = sqrt(sd[j] / (n - 1));

   sample_free(&sm);
}

int main(int argc, char* argv[])
{
   int i;
   char *infilename = NULL;
   char *outfilename = NULL;
   FILE* infile;
   FILE* outfile;
   double *mean, *sd;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-i"))
      {
	 i++;
	 infilename = argv[i];
      }
      else if(strcmp2(argv[i], "-o"))
      {
	 i++;
	 outfilename = argv[i];
      }
   }


   
   


}

