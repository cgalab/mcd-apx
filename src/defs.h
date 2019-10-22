#ifndef MIN_CONVEX_DECOMP_H
#define MIN_CONVEX_DECOMP_H


#ifdef BOOL_DEFINED
typedef bool boolean;
#else
#define false 0
#define true  (!false)
typedef unsigned char  boolean;
#endif

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

#endif
