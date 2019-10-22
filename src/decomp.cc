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


void AddToConvexChain(int *convex, int *num_convex, node *nodes, int start,
                      int end, boolean ccw)
{
   if (ccw) {
      while (start != end) {
         convex[*num_convex] = start;
         ++(*num_convex);
         start = nodes[start].next;
      }
   }
   else {
      while (start != end) {
         convex[*num_convex] = start;
         ++(*num_convex);
         start = nodes[start].prev;
      }
   }

   convex[*num_convex] = start;
   ++(*num_convex);

   return;
}


/*                                                                           */
/* determine the CCW-most vertex that is right of the line  j1-->j2.         */
/*                                                                           */
void GetCCWmostVertexInHalfplane(int j1, int j2, int i_start, int *i2, 
                                 node *nodes, pnt *pnts)
{
   while (CW(&(pnts[j1]), &(pnts[j2]), &(pnts[*i2]))  &&
          (i_start != *i2)) {
      /*                                                                     */
      /* search in CCW direction from *i2.                                   */
      /*                                                                     */
      *i2 = nodes[*i2].next;
   }
   if (!CW(&(pnts[j1]), &(pnts[j2]), &(pnts[*i2])))  
      *i2 = nodes[*i2].prev;

   return;
}


/*                                                                           */
/* determine the CW-most vertex that is right of the line  i1-->i2.          */
/*                                                                           */
void GetCWmostVertexInHalfplane(int i1, int i2, int *i3, node *nodes, 
                                pnt *pnts)
{
   if (CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*i3]))) {
      do {
         /*                                                                  */
         /* search in CW direction from *i3.                                 */
         /*                                                                  */
         *i3  = nodes[*i3].prev;
      } while (CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*i3])));
      *i3  = nodes[*i3].next;
   }
   else { 
      while (!CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*i3]))) {
         *i3  = nodes[*i3].next;
      }
   }

   return;
}


/*                                                                           */
/* determine the nodes  o3  CCW to  o4  that are right of the line  i1-->i2. */
/*                                                                           */
void GetVerticesInHalfplane(int i1, int i2, int *o3, int *o4, node *nodes, 
                            pnt *pnts)
{
   if (CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*o3]))) {
      do {
         /*                                                                  */
         /* search in CW direction from o3.                                  */
         /*                                                                  */
         *o3  = nodes[*o3].prev;
      } while (CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*o3])));
      *o3 = nodes[*o3].next;
   }

   while (!CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*o3]))) {
      *o3 = nodes[*o3].next;
   }

   *o4  = *o3;

   do {
      *o4  = nodes[*o4].next;
   } while (CW(&(pnts[i1]), &(pnts[i2]), &(pnts[*o4])));
   *o4 = nodes[*o4].prev;

   return;
}



int DetermineLowerBound(pnt *pnts, int num_pnts, loop *layers, 
                        int num_layers, node *nodes)
{
   int type_a = 0, type_b = 0;
   int start, curr, next, prev, k, i;
   int outer;

   if (num_layers <= 0)     throw CH_MESSED_UP;
   if (layers[0].num <= 2)  throw CH_MESSED_UP;

   if (num_layers == 1)     return 1;
   if (layers[1].num == 1)  return 3;
   if (layers[1].num == 2)  return 4;

   start = layers[1].nde;
   curr  = start;
   prev  = nodes[curr].prev;
   outer = layers[0].nde;
   do {
      next = nodes[curr].next;
      /*                                                                     */
      /* check whether there exists a vertex on outer CH that is within the  */
      /* cone defined by  curr.                                              */
      /*                                                                     */
      GetCWmostVertexInHalfplane(curr, next, &outer, nodes, pnts);
      if (CW(&(pnts[prev]), &(pnts[curr]), 
             &(pnts[outer]))) {
         ++type_a;
      }
      else {
         ++type_b;
      }
      prev = curr;
      curr = next;
   } while (curr != start);

   k = num_pnts - layers[0].num;
   i = (k + type_a) / 2;
   if ((2 * i) != (k + type_a)) ++i;

   return (1 + i + type_b);
}


void HandleOnePoint(FILE *output, pnt *pnts, loop *layers, node *nodes,
                    int *num_cvx_areas, int L0, int *convex, boolean obj)
{
   int head, next, prev;
   int i1, i2, i3, j1;
   int L1;
   int num_convex = 0;

   L1 = L0 + 1;
   i1 = layers[L0].nde;
   j1 = layers[L1].nde;
   head = i1;

   /*                                                                        */
   /* determine CH vertices that are right of  i1-->j1                       */
   /*                                                                        */
   GetVerticesInHalfplane(i1, j1, &head, &i2, nodes, pnts);
   if (head != i1) printf("head != i1!!!!!\n");
   if (collinear(&(pnts[i1]), &(pnts[j1]), &(pnts[i2])))
      /*                                                                     */
      /* the points  i1, j1, i2  are collinear; only two polygons            */
      /*                                                                     */
      i3 = i2;
   else 
      i3 = nodes[i2].next;
   prev = nodes[i1].prev;
   next = nodes[i3].next;

   AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
   AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
   WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
   ++(*num_cvx_areas);

   if (i2 != i3) {
      /*                                                                     */
      /* the points  i1, j1, i2  are not collinear; three polygons           */
      /*                                                                     */
      AddToConvexChain(convex, &num_convex, nodes, i2, i3, true);
      AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
   }
   
   if (prev == i3) {
      if (next != i1)  printf("next != i1\n");
      AppendNode(nodes, i1, j1);
      AppendNode(nodes, j1, i3);
      CloseLoop(nodes, i1, i3);
      layers[L0].nde = j1;
      layers[L0].num = 3;
   }
   else {
      AppendNode(nodes, prev, i1);
      AppendNode(nodes, i1, j1);
      AppendNode(nodes, j1, i3);
      CloseLoop(nodes, next, i3);
      layers[L0].nde = j1;
      layers[L0].num = CountNodes(nodes, j1);
   }

   return;
}


void HandleTwoPoints(FILE *output, pnt *pnts, loop *layers, node *nodes,
                     int *num_cvx_areas, int L0, int *convex, boolean obj)
{
   int i1, i2, i3, i4, j1, j2, prev, next;
   int L1;
   int num_convex = 0;

   L1 = L0 + 1;
   i1 = layers[L0].nde;
   j1 = layers[L1].nde;
   j2 = nodes[j1].next;

   /*                                                                        */
   /* determine CH vertices that are right of  j1-->j2                       */
   /*                                                                        */
   GetVerticesInHalfplane(j1, j2, &i1, &i2, nodes, pnts);

   if (CCW(&(pnts[i1]), &(pnts[j1]), &(pnts[j2]))) 
      /*                                                                     */
      /* the points  i1, j1, j2  are collinear                               */
      /*                                                                     */
      i4 = i1;
   else
      i4 = nodes[i1].prev;
   if (CCW(&(pnts[j1]), &(pnts[j2]), &(pnts[i2]))) 
      /*                                                                     */
      /* the points  i1, j1, j2  are collinear                               */
      /*                                                                     */
      i3 = i2;
   else
      i3 = nodes[i2].next;
   prev = nodes[i4].prev;
   next = nodes[i3].next;
   AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
   AddToConvexChain(convex, &num_convex, nodes, j2, j1, false);
   WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
   ++(*num_cvx_areas);

   if (i2 != i3) {
      AddToConvexChain(convex, &num_convex, nodes, i2, i3, true);
      AddToConvexChain(convex, &num_convex, nodes, j2, j2, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
   }

   if (i1 != i4) {
      AddToConvexChain(convex, &num_convex, nodes, i4, i1, true);
      AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
   }

   if (i3 == i4) { 
      AppendNode(nodes, i4, j1);
      AppendNode(nodes, j1, j2);
      CloseLoop(nodes, i4, j2);
      layers[L0].nde = i4;
      layers[L0].num = 3;
   }
   else {
      AppendNode(nodes, prev, i4);
      AppendNode(nodes, i4, j1);
      AppendNode(nodes, j1, j2);
      AppendNode(nodes, j2, i3);
      CloseLoop(nodes, next, i3);
      layers[L0].nde = i3;
      layers[L0].num = CountNodes(nodes, i3);
   }

   return;
}


void GeneralCase(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                 pnt *vtx, int *num_cvx_areas, int *convex, boolean obj)
{
   int i1, i2, i3, i4, j1, j2, j3, j0, prev, next, jj;
   int num_convex = 0;

   // test loops
   //i1 = CountNodes(nodes, layers[0].nde);
   //i1 = CountNodes(nodes, layers[1].nde);
   //printf("\n\nafter safety count\n");

   i1 = layers[0].nde;
   j1 = layers[1].nde;
   j2 = nodes[j1].next;
   //PrintNode(nodes, i1, "i1");
   //PrintNode(nodes, j1, "j1");
   //PrintNode(nodes, j2, "j2");

   /*                                                                        */
   /* determine CH vertices that are right of  j1-->j2                       */
   /*                                                                        */
   j3 = nodes[j2].next;
   while (collinear(&(pnts[j1]), &(pnts[j2]), &(pnts[j3]))  &&  (j3 != j1)) 
      j3 = nodes[j3].next;
   if (j3 == j1) {
      /*                                                                     */
      /* all vertices of this (innermost) layer are collinear                */
      /*                                                                     */
      HandleDegenerateLoop(output, pnts, layers, nodes, num_cvx_areas, 0, 
                           convex, obj);
      /*                                                                     */
      /* mark those vertices of the inner CH that end up on the outer CH     */
      /*                                                                     */
      SetConvexHullFlags(vtx, nodes, j1, j2);
      return;
   }
   j2 = nodes[j3].prev;
   GetVerticesInHalfplane(j1, j2, &i1, &i2, nodes, pnts);

   //printf("after GetVerticesInHalfplane()\n");
   //PrintNode(nodes, i2, "i2");

   if (collinear(&(pnts[i1]), &(pnts[j1]), &(pnts[j2]))) 
      /*                                                                     */
      /* the points  i1, j1, j2  are collinear                               */
      /*                                                                     */
      i4 = i1;
   else
      i4 = nodes[i1].prev;
   if (collinear(&(pnts[j1]), &(pnts[j2]), &(pnts[i2]))) 
      /*                                                                     */
      /* the points  j1, j2, j3  are collinear                               */
      /*                                                                     */
      i3 = i2;
   else
      i3 = nodes[i2].next;

   //PrintNode(nodes, i3, "i3");
   //PrintNode(nodes, i4, "i4");

   prev = nodes[j1].prev;
   next = nodes[j2].next;

   /*                                                                        */
   /* output the convex area right of  j1-->j2                               */
   /*                                                                        */
   AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
   AddToConvexChain(convex, &num_convex, nodes, j2, j1, false);
   WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
   ++(*num_cvx_areas);
   //printf("after WriteConvexChain()\n");

   if (i2 != i3) {
      /*                                                                     */
      /* determine all nodes of inner loop right of  jj-->i3                 */
      /*                                                                     */
      jj = j2;
      j3 = nodes[j2].next;
      while (CW(&(pnts[jj]),&(pnts[i3]),&(pnts[j3]))) {
         jj = j3;
         j3 = nodes[j3].next;
      }
      j3 = nodes[j3].prev;
      //printf("after tangent while I\n");
   }
   else {
      j3 = j2;
   }
   //PrintNode(nodes, i3, "i3");   

   if (i1 != i4) {
      /*                                                                     */
      /* determine all nodes of inner loop left of  jj-->i4                  */
      /*                                                                     */
      jj = j1;
      j0 = nodes[j1].prev;
      while (CCW(&(pnts[jj]),&(pnts[i4]),&(pnts[j0]))) {
         jj = j0;
         j0 = nodes[j0].prev;
      }
      j0 = nodes[j0].next;
      //printf("after tangent while II\n");
   }
   else {
      j0 = j1;
   }

   //PrintNode(nodes, j0, "j0");
   //PrintNode(nodes, j3, "j3");   

   /*                                                                        */
   /* assemble new outer and (rest of) inner loop                            */
   /*                                                                        */
   prev = nodes[j0].prev;
   next = nodes[j3].next;
   //PrintNode(nodes, prev, "pr");
   //PrintNode(nodes, next, "ne");   
   CutLoopBeforeAfterNode(nodes, j0, j3);
   //printf("after cut\n");
   SpliceLoops(nodes, prev, next);
   //PrintNode(nodes, prev, "pr");
   //PrintNode(nodes, next, "ne");   
   //PrintNode(nodes, i3, "i3");   
   SpliceLoops(nodes, j3, i3);
   //PrintNode(nodes, j3, "j3");
   //PrintNode(nodes, i3, "i3");   
   SpliceLoops(nodes, i4, j0);
   //PrintNode(nodes, i4, "i4");
   //PrintNode(nodes, j0, "j0");   
   layers[1].nde = next;
   layers[1].num = CountNodes(nodes, next);
   //printf("after CountNodes 1()\n");
   //PrintNode(nodes, i2, "i2");
   //PrintNode(nodes, i3, "i3");
   //PrintNode(nodes, j1, "j1");
   //PrintNode(nodes, j2, "j2");
   layers[0].nde = j1;
   layers[0].num = CountNodes(nodes, j1);
   //printf("after CountNodes 0()\n");
   /*                                                                        */
   /* mark those vertices of the inner CH that end up on the outer CH        */
   /*                                                                        */
   SetConvexHullFlags(vtx, nodes, j0, j3);
   //printf("after SetConvexHullFlags()\n");
   j3 = nodes[j2].next;
   while (j2 != i3) {
      /*                                                                     */
      /* output triangles (partially) right of  j2-->i3                      */
      /*                                                                     */
      AddToConvexChain(convex, &num_convex, nodes, i2, i2, true);
      AddToConvexChain(convex, &num_convex, nodes, j3, j2, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
      j2 = j3;
      j3 = nodes[j3].next;
   }

   j0 = nodes[j1].prev;
   while (j1 != i4) {
      /*                                                                     */
      /* output triangles (partially) left of  j1-->i4                       */
      /*                                                                     */
      AddToConvexChain(convex, &num_convex, nodes, i1, i1, true);
      AddToConvexChain(convex, &num_convex, nodes, j1, j0, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
      j1 = j0;
      j0 = nodes[j0].prev;
   }

   return;
}


void ComputeApproxDecomp(FILE *output, pnt *pnts, int num_pnts,
                         loop *layers, node *nodes,
                         boolean randomized, boolean obj)
{
   pnt *vtx = (pnt*) malloc(MAX * sizeof(pnt));
   int *ch_vtx = (int*) malloc(MAX * sizeof(int));
   int *convex = (int*) malloc(MAX * sizeof(int));
   int i, num_vtx, num_ch_vtx;
   boolean inner_pnts_left = true;
   int num_CH_computations = 0;
   int lower_bound = 0;
   int num_cvx_areas = 0;

   /*                                                                        */
   /* initialize the onion nodes                                             */
   /*                                                                        */
   for (i = 0;  i < num_pnts;  ++i)  {
      nodes[i].vtx  = NIL;
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

   /*                                                                        */
   /* compute and store outermost CH                                         */
   /*                                                                        */
   ConvexHull(vtx, num_vtx, ch_vtx, &num_ch_vtx);
   /*
   //printf("round %d: %d CH vertices out of %d input points\n", 
          num_CH_computations, num_ch_vtx, num_vtx);
   */
   StoreAsOnionLayer(ch_vtx, num_ch_vtx, layers, 0, vtx, nodes, randomized);
   ++num_CH_computations;
   inner_pnts_left = ExtractInnerPoints(vtx, &num_vtx, ch_vtx, num_ch_vtx);

   while (inner_pnts_left) {
      ConvexHull(vtx, num_vtx, ch_vtx, &num_ch_vtx);
      /*
      //printf("round %d: %d CH vertices out of %d input points\n", 
             num_CH_computations, num_ch_vtx, num_vtx);
      */
      StoreAsOnionLayer(ch_vtx, num_ch_vtx, layers, 1, vtx, nodes, randomized);
      ++num_CH_computations;
      if (num_CH_computations == 2) 
         lower_bound = DetermineLowerBound(pnts, num_pnts, layers, 2, nodes);
      if (num_ch_vtx == 1) {
         /*                                                                  */
         /* the innermost layer contains just one point                      */
         /*                                                                  */
         inner_pnts_left = false;
         HandleOnePoint(output, pnts, layers, nodes, &num_cvx_areas, 0, 
                        convex, obj);
      }
      else if (num_ch_vtx == 2) {
         /*                                                                  */
         /* the innermost layer contains just two points                     */
         /*                                                                  */
         inner_pnts_left = false;
         HandleTwoPoints(output, pnts, layers, nodes, &num_cvx_areas, 0, 
                         convex, obj);
      }
      else {
         GeneralCase(output, pnts, layers, nodes, vtx, &num_cvx_areas, convex,
                     obj);

         inner_pnts_left = UpdateInnerPoints(vtx, &num_vtx);
         //if (num_CH_computations == 99) inner_pnts_left = false;
      }
   }

   //ConvexHull(vtx, num_vtx, ch_vtx, &num_ch_vtx);
   //StoreAsOnionLayer(ch_vtx, num_ch_vtx, layers, 1, vtx, nodes, randomized);

   WriteLayers(output, pnts, layers, 1, nodes, obj);
   ++num_cvx_areas;

   if (num_CH_computations == 1) lower_bound = 1;
   printf("#(CH computs): %d\n", num_CH_computations);
   printf("lower bound:   %d\n", lower_bound);
   printf("num_cvx_areas: %d\n", num_cvx_areas);
   printf("apx ratio:     %5.3f\n", ((double) num_cvx_areas) / 
          ((double) lower_bound));
   
   FreeHulls();
   free(vtx);
   free(ch_vtx);
   free(convex);

   return;
}



void HandleDegenerateLoop(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                          int *num_cvx_areas, int L0, int *convex, 
                          boolean obj)
{
   int i1, i2, i3, i4, j1, j2, j3, prev, next;
   int num_convex = 0, L1;

   //printf("in HandleOnionAnnulus(): onion = %d !!!!!!!!!!!!!\n", L0);
   L1 = L0 + 1;
   i1 = layers[L0].nde;
   j1 = layers[L1].nde;
   j2 = nodes[j1].next;
   //PrintNode(nodes, i1, "i1");
   //PrintNode(nodes, j1, "j1");
   //PrintNode(nodes, j2, "j2");

   /*
   if (L0 == 1) {
      i2 = i1;
      printf("outer loop\n");
      do {
         printf("%2d: (%d, %d)\n", i2, (int) pnts[i2].x, (int) pnts[i2].y);
         i2 = nodes[i2].next;
      } while (i2 != i1);
      printf("inner loop\n");
      i2 = j1;
      do {
         printf("%2d: (%d, %d)\n", i2, (int) pnts[i2].x, (int) pnts[i2].y);
         i2 = nodes[i2].next;
      } while (i2 != j1);
   }
   */

   j3 = nodes[j2].next;
   while (collinear(&(pnts[j1]), &(pnts[j2]), &(pnts[j3]))  &&  (j3 > j2))
      j3 = nodes[j3].next;
   j2 = nodes[j3].prev;

   /*                                                                        */
   /* determine CH vertices that are right of  j1-->j2                       */
   /*                                                                        */
   GetVerticesInHalfplane(j1, j2, &i1, &i2, nodes, pnts);

   if (CCW(&(pnts[i1]), &(pnts[j1]), &(pnts[j2]))) 
      /*                                                                     */
      /* the points  i1, j1, j2  are collinear                               */
      /*                                                                     */
      i4 = i1;
   else
      i4 = nodes[i1].prev;
   if (CCW(&(pnts[j1]), &(pnts[j2]), &(pnts[i2]))) 
      /*                                                                     */
      /* the points  i1, j1, j2  are collinear                               */
      /*                                                                     */
      i3 = i2;
   else
      i3 = nodes[i2].next;
   prev = nodes[i4].prev;
   next = nodes[i3].next;
   AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
   AddToConvexChain(convex, &num_convex, nodes, j2, j1, false);
   WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
  ++(*num_cvx_areas);

   if (i2 != i3) {
      AddToConvexChain(convex, &num_convex, nodes, i2, i3, true);
      AddToConvexChain(convex, &num_convex, nodes, j2, j2, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
   }

   if (i1 != i4) {
      AddToConvexChain(convex, &num_convex, nodes, i4, i1, true);
      AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
   }

   if (i3 == i4) { 
      AppendNode(nodes, i4, j1);
      CloseLoop(nodes, i4, j2);
      layers[L0].nde = i4;
      layers[L0].num = CountNodes(nodes, i4);
      layers[L1].nde = i4;
      layers[L1].num = CountNodes(nodes, i4);
   }
   else {
      AppendNode(nodes, prev, i4);
      AppendNode(nodes, i4, j1);
      AppendNode(nodes, j2, i3);
      CloseLoop(nodes, next, i3);
      layers[L0].nde = i4;
      layers[L0].num = CountNodes(nodes, i4);
      layers[L1].nde = i3;
      layers[L1].num = CountNodes(nodes, i3);
   }

   return;

}



void HandleOnionAnnulus(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                        int *num_cvx_areas, int L0, int *convex, boolean obj)
{
   int i1, i2, i3, j1, j2, j3;
   int i_start, j_start;
   int num_convex = 0;

   //printf("in HandleOnionAnnulus(): onion = %d !!!!!!!!!!!!!\n", L0);
   i1 = layers[L0].nde;
   j1 = layers[L0+1].nde;
   j2 = nodes[j1].next;
   //PrintNode(nodes, i1, "i1");
   //PrintNode(nodes, j1, "j1");
   //PrintNode(nodes, j2, "j2");

   /*
   if (L0 == 1) {
      i2 = i1;
      printf("outer loop\n");
      do {
         printf("%2d: (%d, %d)\n", i2, (int) pnts[i2].x, (int) pnts[i2].y);
         i2 = nodes[i2].next;
      } while (i2 != i1);
      printf("inner loop\n");
      i2 = j1;
      do {
         printf("%2d: (%d, %d)\n", i2, (int) pnts[i2].x, (int) pnts[i2].y);
         i2 = nodes[i2].next;
      } while (i2 != j1);
   }
   */

   /*                                                                        */
   /* determine CH vertices that are right of  j1-->j2                       */
   /*                                                                        */
   j3 = nodes[j2].next;
   while (collinear(&(pnts[j1]), &(pnts[j2]), &(pnts[j3]))  &&  (j3 > j2))
      j3 = nodes[j3].next;
   if (j3 < j2) {
      /*                                                                     */
      /* all vertices of this (innermost) layer are collinear                */
      /*                                                                     */
      HandleDegenerateLoop(output, pnts, layers, nodes, num_cvx_areas, L0, 
                           convex, obj);
      return;
   }
   j2 = nodes[j3].prev;
   i_start = NIL;
   j_start = j1;
   GetVerticesInHalfplane(j1, j2, &i1, &i2, nodes, pnts);

   //printf("after GetVerticesInHalfplane()\n");
   //PrintNode(nodes, i1, "i1");
   //PrintNode(nodes, i2, "i2");
   if ((i1 != i2)  &&  (i_start == NIL)) i_start = i1;      

   /*                                                                        */
   /* output the convex area right of  j1-->j2                               */
   /*                                                                        */
   AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
   AddToConvexChain(convex, &num_convex, nodes, j2, j1, false);
   WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
   ++(*num_cvx_areas);
   //printf("after WriteConvexChain()\n");

   j1 = j2;
   j2 = nodes[j2].next;
   //printf("advance j\n");
   //PrintNode(nodes, j1, "j1");
   //PrintNode(nodes, j2, "j2");
   if (!CW(&(pnts[j1]), &(pnts[j2]), &(pnts[i2]))) {
      i3 = nodes[i2].next;
      AddToConvexChain(convex, &num_convex, nodes, i2, i3, true);
      AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      //printf("after WriteConvexChain()\n");
      ++(*num_cvx_areas);
      i2 = i3;
      if (i_start == NIL) i_start = i1;      
   }

   while (j1 != j_start) {
      //printf("************* looping\n");
      i1 = i2;
      //PrintNode(nodes, i1, "i1");   

      /*                                                                     */
      /* determine the next CCW-most  vertex that is right of  j1-->j2       */
      /*                                                                     */
      j3 = nodes[j2].next;
      while (collinear(&(pnts[j1]), &(pnts[j2]), &(pnts[j3]))) 
         j3 = nodes[j3].next;
      j2 = nodes[j3].prev;
      GetCCWmostVertexInHalfplane(j1, j2, i_start, &i2, nodes, pnts);
      //PrintNode(nodes, i2, "i2");   
      if ((i1 != i2)  &&  (i_start == NIL)) i_start = i1;      

      /*                                                                     */
      /* output the convex area right of  j1-->j2                            */
      /*                                                                     */
      AddToConvexChain(convex, &num_convex, nodes, i1, i2, true);
      AddToConvexChain(convex, &num_convex, nodes, j2, j1, false);
      WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
      ++(*num_cvx_areas);
      //printf("after WriteConvexChain()\n");
            
      j1 = j2;
      j2 = nodes[j2].next;
      //printf("advance j\n");
      //PrintNode(nodes, j1, "j1");
      //PrintNode(nodes, j2, "j2");
      if (!CW(&(pnts[j1]), &(pnts[j2]), &(pnts[i2]))  &&  (i2 != i_start)) {
         i3 = nodes[i2].next;
         AddToConvexChain(convex, &num_convex, nodes, i2, i3, true);
         AddToConvexChain(convex, &num_convex, nodes, j1, j1, false);
         WriteConvexChain(output, pnts, nodes, convex, &num_convex, obj);
         //printf("after WriteConvexChain()\n");
         ++(*num_cvx_areas);
         i2 = i3;
         if (i_start == NIL) i_start = i1;      
      }
   }


   return;
}


void ComputeApproxDecompOnion(FILE *output, pnt *pnts, int num_pnts,
                              loop *layers, int num_layers, node *nodes,
                              int lower_bound, boolean obj)
{
   int L0, L1, id;
   int num_cvx_areas = 0;
   int *convex = (int*) malloc(MAX * sizeof(int));

   --num_layers;
   id = num_layers;

   for (L0 = 0;  L0 < num_layers - 0 ;  ++L0) {    
      //printf("num_layers = %d: layer = %d\n", num_layers, L0);

      L1 = L0 + 1;
      if (layers[L1].num == 1) {
         HandleOnePoint(output, pnts, layers, nodes, &num_cvx_areas, L0, 
                        convex, obj);
         id = L0;
      }
      else if (layers[L1].num == 2) {
         HandleTwoPoints(output, pnts, layers, nodes, &num_cvx_areas, L0,
                         convex, obj);
         id = L0;
      }
      else {
         HandleOnionAnnulus(output, pnts, layers, nodes, &num_cvx_areas, L0, 
                            convex, obj);
      }
   }
   /*
   for (L0 = num_layers - 0;  L0 < num_layers;  ++L0) {    
   printf("num_layers = %d: layer = %d  output only\n", num_layers, L0);
   WriteOneLayer(output, pnts, layers, L0, nodes, obj);
   }
   */

   WriteOneLayer(output, pnts, layers, id, nodes, obj);
   ++num_cvx_areas;
      
   printf("lower bound:   %d\n", lower_bound);
   printf("num_cvx_areas: %d\n", num_cvx_areas);
   printf("apx ratio:     %5.3f\n", ((double) num_cvx_areas) /
          ((double) lower_bound));

   free(convex);

   return;
}
