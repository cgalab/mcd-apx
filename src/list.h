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

#ifndef MIN_CONVEX_DECOMP_LIST_H
#define MIN_CONVEX_DECOMP_LIST_H


inline void CutLoopBeforeAfterNode(node* nodes, int head, int tail)
{
   int prev, next;

   prev = nodes[head].prev;
   next = nodes[tail].next;
 
   nodes[prev].next = NIL;
   nodes[next].prev = NIL;

   nodes[head].prev = NIL;
   nodes[tail].next = NIL;
}


inline void CloseLoop(node* nodes, int head, int tail)
{
   nodes[head].prev = tail;
   nodes[tail].next = head;
}


inline void SpliceLoops(node* nodes, int tail, int head)
{
   nodes[tail].next = head;
   nodes[head].prev = tail;
}


inline void AppendNode(node* nodes, int tail, int nde)
{
   nodes[tail].next = nde;
   nodes[nde].prev  = tail;
}


inline int CountNodes(node* nodes, int head)
{
   int num = 0, next;
   
   next = head;

   do {
      ++num;
      next = nodes[next].next;
      //      printf("head = %d, next = %d\n", head, next);
   } while (next != head);

   return num;
}


inline void PrintNode(node* nodes, int nde, char const *s)
{
   printf("%s\n", s);
   printf("nodes[%d].prev = %d\n", nde, nodes[nde].prev);
   printf("nodes[%d].next = %d\n", nde, nodes[nde].next);

   return;
}


#endif
