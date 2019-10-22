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
#include "numerics.h"


void ConvexHull(vertex *vtx, int num_vtx, int *ch_vtx, int *num_ch_vtx)
{
   int *uhull = (int*) malloc(MAX * sizeof(int));
   int *lhull = (int*) malloc(MAX * sizeof(int));
   int i, k, m, uk, lk;

   /* lower hull */
   k = 0;
   for (i = 0;  i < num_vtx;  ++i) {
      while (k >= 2 && !CCW(&(vtx[lhull[k-2]]), &(vtx[lhull[k-1]]), &(vtx[i]))) --k;
      lhull[k++] = i;
   }
   lk = k;

   /* upper hull */
   k = 0;
   for (i = 0;  i < num_vtx;  ++i) {
      while (k >= 2 && !CW(&(vtx[uhull[k-2]]), &(vtx[uhull[k-1]]), &(vtx[i]))) --k;
      uhull[k++] = i;
   }
   uk = k - 2;

   *num_ch_vtx = lk + uk;

   for (i = 0;  i < lk;  ++i) {
      ch_vtx[i] = lhull[i];
      //printf("ch_vtx[%d] = %d = (%f %f)\n", i, ch_vtx[i], vtx[ch_vtx[i]].x, vtx[ch_vtx[i]].y);
   }
   k = lk;
   for (i = uk;  i > 0;  --i) {
      ch_vtx[k++] = uhull[i];
      //printf("ch_vtx[%d] = %d = (%f %f)\n", k-1, ch_vtx[k-1], vtx[ch_vtx[k-1]].x, vtx[ch_vtx[k-1]].y);
   }

   /*
   for (i = 0;  i < *num_ch_vtx;  ++i)
      printf("ch_vtx[%d] = %d = (%f %f)\n", i, ch_vtx[i], 
             vtx[ch_vtx[i]].x, vtx[ch_vtx[i]].y);
   */

   free(lhull);
   free(uhull);

   return;
}


void GetInnerPnts(vertex *vtx, int num_vtx, int *ch_vtx, int num_ch_vtx,
                  int *in_vtx, int *num_in_vtx)
{
   int *pnts  = (int*) malloc(MAX * sizeof(int));
   int i, j = 0;

   for (i = 0; i < num_vtx;  ++i) {
      pnts[i] = i;
   }
   
   for (i = 0;  i < num_ch_vtx;  ++i) pnts[ch_vtx[i]] = -1;
   
   for (i = 0;  i < num_vtx;  ++i) {
      if (pnts[i] >= 0) {
         in_vtx[j++] = pnts[i];
      }
   }
   *num_in_vtx = j;

   /*
   for (i = 0;  i < *num_in_vtx;  ++i)
      printf("in_vtx[%d] = %d = (%f %f)\n", i, in_vtx[i], 
             vtx[in_vtx[i]].x, vtx[in_vtx[i]].y);
   */
   (void) vtx;
   
   free(pnts);

   return;
}
