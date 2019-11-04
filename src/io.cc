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

#include <math.h>

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



/*
 * INPUT FOMRAT:
 * # as comments any line starting with a # is ignored
 * vertices are given as vertexidx x-coordinate y-coordinate
 * -- extension: if only two values are given then x/y coordinates
 *    are assumed and the index is simple counted line by line
 * */
void ReadInput(FILE *input, pnt *pnts, int *num_pnts) {
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	int cnt = 0;

	if (input == NULL) {throw EOF_ENCOUNTERED;}

	while ((read = getline(&line, &len, input)) != -1) {
		if(read > 0) {
			if(line[0] == '#' || line[0] == '\n') /* comment, NL */ {continue;}
			else {
				double id = NIL;
				if (2 == sscanf(line, "%lf %lf %lf\n", &id, &pnts[cnt].x, &pnts[cnt].y)) {
					pnts[cnt].y = pnts[cnt].x;
					pnts[cnt].x = id;
					pnts[cnt].id = cnt;
				} else {
					pnts[cnt].id = (int)id;
				}
				pnts[cnt].in = true;
				++cnt;
			}
		}
	}

	*num_pnts = cnt;

//   int i = 0;
//
//   if (EOF == fscanf(input, "%*[^\n]\n")) throw EOF_ENCOUNTERED;
//   if (EOF == fscanf(input, "%*[^\n]\n")) throw EOF_ENCOUNTERED;
//
//   if (EOF == fscanf(input, "%d %lf %lf",
//                     &(pnts[0].id), &(pnts[0].x), &(pnts[0].y)))
//      throw EOF_ENCOUNTERED;
//
//   if (EOF == fscanf(input, "%d %lf %lf",
//                     &(pnts[1].id), &(pnts[1].x), &(pnts[1].y)))
//      throw EOF_ENCOUNTERED;
//
//   if (EOF == fscanf(input, "%d %lf %lf",
//                     &(pnts[2].id), &(pnts[2].x), &(pnts[2].y)))
//      throw EOF_ENCOUNTERED;
//
//   pnts[0].in = pnts[1].in = pnts[2].in = true;
//
//   i = 3;
//   while (EOF != fscanf(input, "%d %lf %lf",
//                        &(pnts[i].id), &(pnts[i].x), &(pnts[i].y))) {
//      pnts[i].in = true;
//      ++i;
//   }
//
//   *num_pnts = i;

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



boolean WriteConvexChain(FILE *output, pnt *pnts, node *nodes, int *convex, 
                         int *num_convex, boolean obj)
{
   int i, j = 0, k = 1;

   for (i = 1;  i < *num_convex; ++i) {
      if (convex[i] != convex[j]) {
         convex[k] = convex[i];
         ++j;
         ++k;
      }
   }
   *num_convex = k;

   if (*num_convex >= 3) {
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
      return true;
   }
   else {
      *num_convex = 0;
      return false;
   }
}



void WriteObjVertices(FILE *output, pnt *pnts, int num_pnts)
{
   int i;

   for (i = 0; i < num_pnts;  ++i) {
	   if(floor(pnts[i].x) == pnts[i].x && floor(pnts[i].y) == pnts[i].y) {
		   fprintf(output, "v %d %d 0\n", (int) pnts[i].x, (int) pnts[i].y);
	   } else {
		   fprintf(output, "v %lf %lf 0\n", pnts[i].x, pnts[i].y);
	   }
   }

   return;
}
