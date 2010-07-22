#include "cd.h"
#include "util.h"

/*
 * Converts a matrix encoded as chars (one byte per entry), scales it, and
 * saves it as a matrix of doubles (8 bytes per entry).
 */
int scale(gmatrix *g, char* filename)
{
   int i, j;
   FILE *fout = NULL, *fin = NULL;
   dtype *tmp;
   double *tmp2;
   double delta;

   MALLOCTEST(g->mean, sizeof(double) * (g->p + 1))
   MALLOCTEST(g->sd, sizeof(double) * (g->p + 1))
   g->mean[0] = 0;
   g->sd[0] = 1;

   MALLOCTEST(tmp, sizeof(dtype) * g->n)
   MALLOCTEST(tmp2, sizeof(double*) * g->n)

   FOPENTEST(fin, g->filename, "rb")
   FOPENTEST(fout, filename, "w")

   /* read y but do not scale it */
   FREADTEST(tmp, sizeof(dtype), g->n, fin)

   for(i = 0 ; i < g->n ; i++)
      tmp2[i] = (double)tmp[i];

   FWRITETEST(tmp2, sizeof(double), g->n, fout)

   /* read the data and scale each variable */
   for(j = 1 ; j < g->p + 1 ; j++)
   {
      printf("%d of %d\r", j, g->p);
      FREADTEST(tmp, sizeof(dtype), g->n, fin)

      g->mean[j] = g->sd[j] = 0;
      for(i = 0 ; i < g->n ; i++)
      {
	 tmp2[i] = (double)tmp[i];

	 delta = tmp2[i] - g->mean[j];
	 g->mean[j] += delta / (i + 1);
	 g->sd[j] += delta * (tmp2[i] - g->mean[j]);
      }

      g->sd[j] = sqrt(g->sd[j] / (g->n - 1));

      for(i = 0 ; i < g->n ; i++)
      {
	 tmp2[i] = tmp2[i] - g->mean[j];
	 if(g->sd[j] > SDTHRESH)
	    tmp2[i] /= g->sd[j];
      }
        
      FWRITETEST(tmp2, sizeof(double), g->n, fout)
   }
   printf("\n");

   fclose(fin);
   fclose(fout);

   free(tmp);
   free(tmp2);
      
   return SUCCESS;
}

int main(int argc, char* argv[])
{
   int i, n = 0, p = 0;
   char *filename_in = NULL;
   char *filename_out = NULL;
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
   }

   if(filename_in == NULL || filename_out == NULL || n == 0 || p == 0)
   {
      printf("scale: -fin <filein> -fout <fileout> -n #n -p #p\n");
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, filename_in, n, p))
      return EXIT_FAILURE;

   scale(&g, filename_out);

   writevectorf("mean.csv", g.mean, p + 1);
   writevectorf("sd.csv", g.sd, p + 1);

   gmatrix_free(&g);

   return EXIT_SUCCESS;
}

