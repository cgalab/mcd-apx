/*****************************************************************************/
/*                                                                           */
/*                        MinConvexDecomp                                    */
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


/*                                                                           */
/* this function parses the command-line arguments                           */
/*                                                                           */
void ArgEval(int argc, char *argv[], rt_options *rt_opt)
{
   int count = 1;
   boolean success = true;

   rt_opt->seed = NIL;
   rt_opt->counter = 1;
   rt_opt->timeout = 0;
   rt_opt->verbose = false;
   rt_opt->index = false;

   /*                                                                        */
   /* parse the command-line arguments                                       */
   /*                                                                        */
   while ((count < argc) && success)    {
      
      if (strcmp(argv[count],"--input") == 0) {
         ++count;
         if ((success = (count < argc))) 
            rt_opt->input_file = argv[count];
      }
      else if (strcmp(argv[count],"--output") == 0) {
         ++count;
         if ((success = (count < argc))) 
            rt_opt->output_file = argv[count];
      }
      else if (strcmp(argv[count],"--seed") == 0) {
         ++count;
         if ((success = (count < argc)))  
            rt_opt->seed = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--counter") == 0) {
         ++count;
         if ((success = (count < argc)))  
            rt_opt->counter = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--timeout") == 0) {
         ++count;
         if ((success = (count < argc)))
            rt_opt->timeout = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--verbose") == 0) {
         rt_opt->verbose = true;
      }
      else if (strcmp(argv[count],"--random") == 0) {
         rt_opt->randomized = true;
      }
      else if (strcmp(argv[count],"--index") == 0) {
         rt_opt->index = true;
      }
      else if (strcmp(argv[count],"--onion") == 0) {
         rt_opt->onion = true;
      }
      else if (strcmp(argv[count],"--obj") == 0) {
         rt_opt->obj = true;
      }
      else {
         throw CL_ARG_ERROR;
      }
      ++count;
   }

   return;
}
