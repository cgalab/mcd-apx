/*****************************************************************************/
/*                                                                           */
/*             Minimum Convex Deomposition via Onion Layers                  */
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
#include <time.h>
#include <cmath>

#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

/*                                                                           */
/* get my header files                                                       */
/*                                                                           */
#include "defs.h"
#include "headers.h"

int p_comp(const void *, const void *);


int main(int argc, char *argv[])
{
   rt_options rt_opt;
   FILE *input, *output;
   pnt *pnts = (pnt*) malloc(MAX * sizeof(pnt));
   loop *layers = (loop*) malloc(MAX * sizeof(loop));
   node *nodes = (node*) malloc(MAX * sizeof(node));
   int num_pnts = 0;
   int num_layers = 0;
   int err = 0;
   int lower_bound;

   try {
      /*                                                                     */
      /* parse command-line arguments                                        */
      /*                                                                     */
      ArgEval(argc, argv, &rt_opt);

      if (rt_opt.randomized) {
         /*                                                                  */
         /* prepare for randomized cutting of convex areas                   */
         /*                                                                  */
         if (rt_opt.seed == NIL) {
            long rseed;
            int need = sizeof(rseed);
#define RND_SOURCE "/dev/urandom"
            int fd = open(RND_SOURCE, O_RDONLY);
            if (fd < 0) {
               fprintf(stderr, "Cannot open " RND_SOURCE ": %s\n", strerror(errno));
               exit(1);
            }
            int cnt = read(fd, &rseed, need);
            if (cnt <= 0) {
               fprintf(stderr, "Error reading from " RND_SOURCE ": %s\n", strerror(errno));
               exit(1);
            }
            if (cnt != need) {
               fprintf(stderr, "Warning: short read from /dev/urandom, randomness not seeded well.\n");
            }
            close(fd);
            rt_opt.seed = rseed;
         }
         printf("random_seed: %ld\n", rt_opt.seed);
         srand48(rt_opt.seed);
      }
         
      /*                                                                     */
      /* read input pnts                                                     */
      /*                                                                     */
      input = OpenFile(rt_opt.input_file.c_str(), "r");
      ReadInput(input, pnts, &num_pnts);
      fclose(input);

      output = OpenFile(rt_opt.output_file.c_str(), "w");

      if (rt_opt.obj) {
         /*                                                                  */
         /* write points for OBJ output                                      */
         /*                                                                  */
         WriteObjVertices(output, pnts, num_pnts);
      }

      /*                                                                     */
      /* sort points in lexicographical order (first x, second y)            */
      /*                                                                     */
      qsort(pnts, num_pnts, sizeof(pnt), p_comp);

      /*
      for (i = 0;  i < num_pnts;  ++i) {
         printf("%d: %f %f\n", pnts[i].id, pnts[i].x, pnts[i].y);
      }
      */

      if (rt_opt.onion) {
         /*                                                                  */
         /* compute onion layers                                             */
         /*                                                                  */
         OnionLayers(pnts, num_pnts, layers, &num_layers, nodes);

         /*                                                                  */
         /* output layers                                                    */
         /*                                                                  */
         /*
         output = OpenFile(rt_opt.output_file, "w");
         WriteLayers(output, pnts, layers, num_layers, nodes);
         fclose(output);
         exit(1);
         */

         /*                                                                  */
         /* compute lower bound                                              */
         /*                                                                  */
         lower_bound = DetermineLowerBound(pnts, num_pnts, layers, num_layers, 
                                           nodes);

         /*                                                                  */
         /* compute approximate minimum decomposition (based on onions)      */
         /*                                                                  */
         ComputeApproxDecompOnion(output, pnts, layers, num_layers,
                                  nodes, lower_bound, rt_opt.obj);
         
      } else if(rt_opt.randomized) {
         /*                                                                  */
         /* compute approximate minimum decomposition (Knauer&Spillner)      */
         /*                                                                  */
         ComputeApproxDecomp(output, pnts, num_pnts, layers, nodes,
                             rt_opt.randomized, rt_opt.obj);
      }
      fclose(output);
   }
   catch (errordef PolyAreaErrorCode) {
      err = 1;
      switch (PolyAreaErrorCode) {
      case SUCCESS: 
         break;
      case CL_ARG_ERROR:
         fprintf(stderr, 
                 "\n*** exception: incorrect command-line arguments ***\n");
         break;
      case MEM_ALLOC_FAILED:
         fprintf(stderr, 
                 "\n*** exception: memory (re-)allocation failed ***\n");
         break;
      case FILE_ACCESS_FAILED: 
         fprintf(stderr, 
                 "\n*** exception: I/O file could not be opened ***\n");
         break;
      case INSUFFICENT_INPUT:
         fprintf(stderr, 
                 "\n*** exception: less than three input vertices ***\n");
         break;
      case EOF_ENCOUNTERED:
         fprintf(stderr, 
                 "\n*** exception: end-of-file encountered during file input ***\n");
         break;
      case INDEX_MISMATCH:
         fprintf(stderr, 
                 "\n*** exception: index mismatch during file input ***\n");
         break;
      case LINK_MISMATCH:
         fprintf(stderr, 
                 "\n*** exception: link mismatch in 2-Opt ***\n");
         break;
      case NUMBER_MISMATCH:
         fprintf(stderr, 
                 "\n*** exception: numbers of CH and inner vertices do not add up ***\n");
         break;
      case LIST_MESSED_UP:
         fprintf(stderr, 
                 "\n*** exception: circular list is messed up ***\n");
         break;
      case CH_MESSED_UP:
         fprintf(stderr, 
                 "\n*** exception: convex hull is messed up ***\n");
         break;
      case ONION_MESSED_UP:
         fprintf(stderr, 
                 "\n*** exception: onion layers are messed up ***\n");
         break;
      case UNKNOWN_ERROR:
         fprintf(stderr, 
                 "\n*** exception: unknown error code ***\n");
         break;
      }
   }

   free(pnts);
   free(nodes);
   free(layers);

   exit(err);
}


int p_comp(const void *a, const void *b)
{
   if      (((pnt*)a)->x < ((pnt*)b)->x)     return -1;
   else if (((pnt*)a)->x > ((pnt*)b)->x)     return  1;
   else  {
      if      (((pnt*)a)->y < ((pnt*)b)->y)  return -1;
      else if (((pnt*)a)->y > ((pnt*)b)->y)  return  1;
      else                                   return  0;
   }
}

