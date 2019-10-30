#ifndef MIN_CONVEX_DECOMP_H
#define MIN_CONVEX_DECOMP_H

#include <assert.h>

#include <list>

using boolean = bool;

typedef struct {
   double x;              /* x coordinate */
   double y;              /* y coordinate */ 
   int id;                /* original id of pnt */
   boolean in;            /* is in interior of CHs */
} pnt;

typedef struct {
   int nde;
   int num;
} loop;


typedef struct {
   int vtx;               /* index in vtx[]/pnts[]  */
   int next;              /* next CCW node          */
   int prev;              /* next CW node           */ 
} node;


typedef struct {
   boolean verbose;
   boolean index;
   int seed;
   int counter;
   int timeout;
   char *input_file;
   char *output_file;
   boolean randomized;
   boolean onion;
   boolean obj;
   int partition;
} rt_options;


typedef enum {
   SUCCESS,
   CL_ARG_ERROR,
   MEM_ALLOC_FAILED,
   FILE_ACCESS_FAILED,
   INSUFFICENT_INPUT,
   EOF_ENCOUNTERED,
   INDEX_MISMATCH,
   LINK_MISMATCH,
   NUMBER_MISMATCH,
   LIST_MESSED_UP,
   CH_MESSED_UP,
   ONION_MESSED_UP,
   UNKNOWN_ERROR
} errordef;

#define MAX  1000001
#define NIL       -1

class Data {
public:
	Data(int k = 1) {
		assert(k>0);
		pnts   = new pnt[MAX/k];
		layers = new loop[MAX/k];
		nodes  = new node[MAX/k];
	}
	~Data() {
		delete[] pnts;
		delete[] nodes;
		delete[] layers;
	}

	int num_pnts    = 0;
	int num_layers  = 0;
	int num_nodes   = 0;
	int lower_bound = 0;

	pnt  *pnts;
	loop *layers;
	node *nodes;
};

#endif
