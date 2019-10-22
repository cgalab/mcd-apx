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
#include <stdlib.h>

/*                                                                           */
/* get my header files                                                       */
/*                                                                           */
#include "defs.h"
#include "headers.h"
#include "list.h"
#include "numerics.h"
#include "random.h"

int *uhull = (int*) malloc(MAX * sizeof(int));
int *lhull = (int*) malloc(MAX * sizeof(int));


void ConvexHull(pnt *vtx, int num_vtx, int *ch_vtx, int *num_ch_vtx)
{
   int i, k, uk, lk;

   if (num_vtx <= 2) {
      ch_vtx[0] = 0;
      if (num_vtx == 2) ch_vtx[1] = 1;
      *num_ch_vtx = num_vtx;
   }
   else {
      
      /* lower hull */
      k = 0;
      for (i = 0;  i < num_vtx;  ++i) {
         while (k >= 2 && 
                !CCW(&(vtx[lhull[k-2]]), &(vtx[lhull[k-1]]), &(vtx[i]))) --k;
         lhull[k++] = i;
      }
      lk = k;
      
      /* upper hull */
      if (lk < num_vtx) {
         k = 0;
         for (i = 0;  i < num_vtx;  ++i) {
            while (k >= 2 && 
                   !CW(&(vtx[uhull[k-2]]), &(vtx[uhull[k-1]]), &(vtx[i]))) --k;
            uhull[k++] = i;
         }
         uk = k - 2;
      }
      else {
         /*                                                                   */
         /* all points are collinear                                          */
         /*                                                                   */
         uk = 0;
      }

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
   }

   return;
}


void StoreAsOnionLayer(int *ch_vtx, int num_ch_vtx, 
                       loop *layers, int layer_id, 
                       pnt *vtx, node *nodes, boolean randomized)
{
   int i, id, prev;
   
   id = vtx[ch_vtx[0]].id;
   layers[layer_id].nde = id;
   layers[layer_id].num = num_ch_vtx;
   nodes[id].vtx  = ch_vtx[0];
   nodes[id].prev = id;
   nodes[id].next = id;
   prev = id;
   for (i = 1;  i < num_ch_vtx;  ++i) {
      id = vtx[ch_vtx[i]].id;
      nodes[id].vtx    = ch_vtx[i];
      nodes[id].prev   = prev;
      nodes[prev].next = id;
      prev = id;
   }
   nodes[prev].next = layers[layer_id].nde;
   nodes[layers[layer_id].nde].prev = prev;

   if (randomized) {
      i = RandomInt(num_ch_vtx);
      layers[layer_id].nde = vtx[ch_vtx[i]].id;
   }

   return;
}


boolean ExtractInnerPoints(pnt *vtx, int *num_vtx, 
                           int *ch_vtx, int num_ch_vtx)
{
   int i, k = 0;

   if (*num_vtx > num_ch_vtx) {
      for (i = 0;  i < num_ch_vtx;  ++i) vtx[ch_vtx[i]].in = false;
      for (i = 0; i < *num_vtx;  ++i) {
         if (vtx[i].in) {
            vtx[k] = vtx[i];
            ++k;
         }
      }
      *num_vtx = k;
   }

   return (k > 0);
}


boolean UpdateInnerPoints(pnt *vtx, int *num_vtx)
{
   int i, k = 0;

   for (i = 0; i < *num_vtx;  ++i) {
      if (vtx[i].in) {
         vtx[k] = vtx[i];
         ++k;
      }
   }
   *num_vtx = k;

   return (k > 0);
}



void SetConvexHullFlags(pnt *vtx, node *nodes, int start, int end)
{
   int nde; 

   nde = start;
   vtx[nodes[nde].vtx].in = false;
   
   while (nde != end) {
      nde = nodes[nde].next;
      vtx[nodes[nde].vtx].in = false;
   }

   return;
}


void FreeHulls(void)
{
   free(lhull);
   free(uhull);

   return;
}


void OnionLayers(pnt *pnts, int num_pnts, loop *layers, int *num_layers, 
                 node *nodes, int *num_nodes)
{
   pnt *vtx = (pnt*) malloc(MAX * sizeof(pnt));
   int *ch_vtx = (int*) malloc(MAX * sizeof(int));
   int i, num_vtx, num_ch_vtx;
   boolean inner_pnts_left = true;

   /*                                                                        */
   /* initialize the onion nodes                                             */
   /*                                                                        */
   for (i = 0;  i < num_pnts;  ++i)  {
      nodes[i].vtx  = i;
      nodes[i].prev = NIL;
      nodes[i].next = NIL;
   }

   /*                                                                        */
   /* make a working copy of the input points                                */
   /*                                                                        */
   num_vtx = num_pnts;
   for (i = 0;  i < num_pnts;  ++i)  {
      vtx[i].x  = pnts[i].x;
      vtx[i].y  = pnts[i].y;
      vtx[i].id = i;
      vtx[i].in = true;
   }

   do {
      /*                                                                     */
      /* compute convex hull of current point set                            */
      /*                                                                     */
      ConvexHull(vtx, num_vtx, ch_vtx, &num_ch_vtx);
      /*
      printf("round %d: %d CH vertices out of %d input points\n", 
             *num_layers, num_ch_vtx, num_vtx);
             */
    
      /*                                                                     */
      /* store current CH as onion layer                                     */
      /*                                                                     */
      StoreAsOnionLayer(ch_vtx, num_ch_vtx, layers, *num_layers, vtx, nodes,
                        false);
      ++(*num_layers);

      /*                                                                     */
      /* extract inner pnts (and discard CH vertices)                        */
      /*                                                                     */
      inner_pnts_left = ExtractInnerPoints(vtx, &num_vtx, ch_vtx, num_ch_vtx);

   } while (inner_pnts_left);
   
   free(vtx);
   free(ch_vtx);
   FreeHulls();

   return;
}

