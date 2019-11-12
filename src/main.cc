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
#include "broker.h"
#include "headers.h"


int main(int argc, char *argv[])
{
	rt_options rt_opt;
	FILE *input = NULL, *output = NULL;

	/*                                                                     */
	/* parse command-line arguments                                        */
	/*                                                                     */
	ArgEval(argc, argv, &rt_opt);

	if (rt_opt.randomized) {
		/*                                                                  */
		/* prepare for randomized cutting of convex areas                   */
		/*                                                                  */
		if (rt_opt.seed == NIL) {
			initRand(&rt_opt);
		}
		printf("random_seed: %ld\n", rt_opt.seed);
		fflush(stdout);
		srand48(rt_opt.seed);
	}

	Broker data;
	data.cfg = &rt_opt;

	/*                                                                     */
	/* read input pnts                                                     */
	/*                                                                     */
	input = OpenFile(rt_opt.input_file.c_str(), "r");
	ReadInput(input, data.pnts, &data.num_pnts);
	fclose(input);

	output = OpenFile(rt_opt.output_file.c_str(), "w");

	if (rt_opt.obj) {
		/*                                                                  */
		/* write points for OBJ output                                      */
		/*                                                                  */
		WriteObjVertices(output, data.pnts, data.num_pnts);
	}

	/*                                                                     */
	/* sort points in lexicographical order (first x, second y)            */
	/*                                                                     */
	qsort(data.pnts, data.num_pnts, sizeof(pnt), p_comp);


	if(rt_opt.partition > 1) {
		/************************************************************/
		/*			      PARTITION APPROACH                      */
		/************************************************************/

		/* split into 'rt_opt.partition' point sets */
		data.partition(rt_opt.partition);

		/* compute onions of each set */
		for(auto& s : data.sets) {
			StartComputation(&s,rt_opt,output);
		}

		/* merge between sets */
		data.merge();

		/* write additional faces from merge*/
		data.writeFacesToFile(output);

	} else if(rt_opt.partition < 1) {
		/************************************************************/
		/*			      RANDOM   APPROACH                      */
		/************************************************************/
		data.runTriangleOnlyApproach();
		data.writeFacesToFile(output);
	} else {
		/************************************************************/
		/*			     CLASSICAL APPROACH                       */
		/************************************************************/
		StartComputation(&data,rt_opt,output);
	}

	/* print stats for wrapper */
	int cvx_faces = 0;
	cvx_faces += data.faces.size();
	for(auto& s : data.sets) {
		cvx_faces += s.num_cvx_areas;
	}
	std::cout << "num_cvx_areas: " << cvx_faces << std::endl;
	std::cout << "num_partitions: " << rt_opt.partition << std::endl;
	std::cout << "num_tri_fips: " << rt_opt.flip_tris << std::endl;
	fflush(stdout);


	FreeHulls();
	fclose(output);

	exit(0);
}


void StartComputation(Data *data, rt_options &rt_opt, FILE *output) {
	if (rt_opt.onion) {
		/*                                                                  */
		/* compute onion layers                                             */
		/*                                                                  */
		OnionLayers(data->pnts, data->num_pnts, data->layers, &data->num_layers, data->nodes);
		/*                                                                  */
		/* output layers                                                    */
		/*                                                                  */

		/*                                                                  */
		/* compute lower bound                                              */
		/*                                                                  */
		data->lower_bound = DetermineLowerBound(data->pnts, data->num_pnts, data->layers, data->num_layers,
				data->nodes);

		if(rt_opt.partition > 1) {
			data->backupOnionZero();
		}

		/*                                                                  */
		/* compute approximate minimum decomposition (based on onions)      */
		/*                                                                  */
		ComputeApproxDecompOnion(data,output, data->pnts, data->num_pnts, data->layers, data->num_layers,
				data->nodes, data->lower_bound, rt_opt.obj, rt_opt);
	}

	if(rt_opt.randomized) {
		/*                                                                  */
		/* compute approximate minimum decomposition (Knauer&Spillner)      */
		/*                                                                  */
		ComputeApproxDecomp(data, output, data->pnts, data->num_pnts, data->layers, data->nodes,
				rt_opt.randomized, rt_opt.obj, rt_opt);
	}

	if(!rt_opt.randomized && !rt_opt.onion) {
		printf("Error: no approach defined!\n");
		exit(0);
	}
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

int p_comp_y(const void *a, const void *b)
{
   if      (((pnt*)a)->y < ((pnt*)b)->y)     return -1;
   else if (((pnt*)a)->y > ((pnt*)b)->y)     return  1;
   else  {
      if      (((pnt*)a)->x < ((pnt*)b)->x)  return -1;
      else if (((pnt*)a)->x > ((pnt*)b)->x)  return  1;
      else                                   return  0;
   }
}

bool p_comp_y(pnt a, pnt b)
{
   if      (a.y < b.y)     return false;
   else if (a.y > b.y)     return true;
   else  {
      if      (a.x < b.x)  return false;
      else if (a.x > b.x)  return true;
      else                 return true;
   }
}

void initRand(rt_options *rt_opt) {
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
	/* fprintf(stderr, "Seeding with %lx\n", rseed); */
	rt_opt->seed = rseed;
}
