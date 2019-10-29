/*****************************************************************************/
/*                                                                           */
/* Written by:  Martin Held                                                  */
/*                                                                           */
/* E-Mail:      held@cs.sbg.ac.at                                            */
/* Fax Mail:    (+43 662) 8044-611                                           */
/* Voice Mail:  (+43 662) 8044-6304                                          */
/* Snail Mail:  Martin Held                                                  */
/*              FB Computerwissenschaften                                    */
/*              Universitaet Salzburg                                        */
/*              A-5020 Salzburg, Austria                                     */
/*                                                                           */
/*****************************************************************************/

/*                                                                           */
/* get standard libraries                                                    */
/*                                                                           */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*                                                                           */
/* get my header files                                                       */
/*                                                                           */
#include "defs.h"
#include "headers.h"


FILE *OpenFile(const char *file_name, const char *access)
{
   FILE *file;
 
   if ((file = fopen(file_name, access)) == NULL) throw FILE_ACCESS_FAILED;

   return file;
}


void ReadInput(FILE *input, pnt *pnts, int *num_pnts)
{
   int i = 0;
   
   if (EOF == fscanf(input, "%*[^\n]\n")) throw EOF_ENCOUNTERED;
   if (EOF == fscanf(input, "%*[^\n]\n")) throw EOF_ENCOUNTERED;
   
   if (EOF == fscanf(input, "%d %lf %lf", 
                     &(pnts[0].id), &(pnts[0].x), &(pnts[0].y)))
      throw EOF_ENCOUNTERED;

   if (EOF == fscanf(input, "%d %lf %lf", 
                     &(pnts[1].id), &(pnts[1].x), &(pnts[1].y)))
      throw EOF_ENCOUNTERED;

   if (EOF == fscanf(input, "%d %lf %lf", 
                     &(pnts[2].id), &(pnts[2].x), &(pnts[2].y)))
      throw EOF_ENCOUNTERED;

   pnts[0].in = pnts[1].in = pnts[2].in = true;

   i = 3;
   while (EOF != fscanf(input, "%d %lf %lf", 
                        &(pnts[i].id), &(pnts[i].x), &(pnts[i].y))) {
      pnts[i].in = true;
      ++i;
   }

   *num_pnts = i;

   for(int j = 0; j<i ; ++j ) {
	   printf("%d %d %lf %lf\n",j,pnts[j].id,pnts[j].x,pnts[j].y);
   }

  return;
}


void WriteOneLayer(FILE *output, pnt *pnts, loop *layers, int layer_id,
                   node *nodes, boolean obj)
{
   int start, next;

   if (obj) {
      fprintf(output, "f");
      if (layers[layer_id].num == 1) {
         start = layers[layer_id].nde;
         fprintf(output, " %d", pnts[start].id + 1);
      }
      else if (layers[layer_id].num == 2) {
         start = layers[layer_id].nde;
         fprintf(output, " %d", pnts[start].id + 1);
         next = nodes[start].next;
         fprintf(output, " %d", pnts[next].id + 1);
      }
      else if (layers[layer_id].num > 2) {
         start = layers[layer_id].nde;
         next  = start;
         do {
            fprintf(output, " %d", pnts[next].id + 1);
            next = nodes[next].next;
         } while (next != start);
      }
      fprintf(output, "\n");
   }
   else {
   //printf("writing layer %d\n", layer_id);
      if (layers[layer_id].num == 1) {
         fprintf(output, "1\n");
         start = layers[layer_id].nde;
         fprintf(output, "%f %f\n", pnts[start].x, pnts[start].y);
      }
      else if (layers[layer_id].num == 2) {
         fprintf(output, "2\n");
         start = layers[layer_id].nde;
         fprintf(output, "%f %f\n", pnts[start].x, pnts[start].y);
         next = nodes[start].next;
         fprintf(output, "%f %f\n", pnts[next].x, pnts[next].y);
      }
      else if (layers[layer_id].num > 2) {
         fprintf(output, "%d\n", layers[layer_id].num+1);
         start = layers[layer_id].nde;
         next  = start;
         do {
            //printf("%d\n", next);
            fprintf(output, "%f %f\n", pnts[next].x, pnts[next].y);
            next = nodes[next].next;
         } while (next != start);
         //printf("%d\n", next);
         fprintf(output, "%f %f\n", pnts[next].x, pnts[next].y);
      }
   }

   return;
}


void WriteOneTriLayer(FILE *output, pnt *pnts, int tri[], node *nodes,
                      boolean obj)
{
   int i;

   if (obj) {
      fprintf(output, "f");
      for (i = 0;  i < 3;  ++i) {
         fprintf(output, " %d", pnts[tri[i]].id + 1);
      }
      fprintf(output, "\n");
   }
   else {
      fprintf(output, "4\n");
      
      for (i = 0;  i < 3;  ++i) {
         fprintf(output, "%f %f\n", pnts[tri[i]].x, pnts[tri[i]].y);
      }
      fprintf(output, "%f %f\n", pnts[tri[0]].x, pnts[tri[0]].y);
   }

   //printf("tri: %d %d %d \n", tri[0], tri[1], tri[2]);

   return;
}



void WriteLayers(FILE *output, pnt *pnts, loop *layers, int num_layers, 
                 node *nodes, boolean obj)
{
   int i;

   for (i = 0;  i < num_layers;  ++i) {
      WriteOneLayer(output, pnts, layers, i, nodes, obj);
   }

   return;
}



void WriteConvexChain(FILE *output, pnt *pnts, node *nodes, int *convex, 
                      int *num_convex, boolean obj)
{
   int i;

   if (obj) {
      fprintf(output, "f");
      for (i = 0;  i < *num_convex;  ++i) {
         fprintf(output, " %d", pnts[convex[i]].id + 1);
      }
      fprintf(output, "\n");
   }
   else {
      fprintf(output, "%d\n", *num_convex + 1);
      
      for (i = 0;  i < *num_convex;  ++i) {
         fprintf(output, "%f %f\n", pnts[convex[i]].x, pnts[convex[i]].y);
      }
      fprintf(output, "%f %f\n", pnts[convex[0]].x, pnts[convex[0]].y);
   }

   *num_convex = 0;

   return;
}



void WriteObjVertices(FILE *output, pnt *pnts, int num_pnts)
{
   int i;

   for (i = 0; i < num_pnts;  ++i) {
      fprintf(output, "v %lf %lf 0.0\n", pnts[i].x, pnts[i].y);
   }

   return;
}
