#ifndef MIN_CONVEX_DECOMP_H
#define MIN_CONVEX_DECOMP_H

#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

#define MAX  1000001
#define NIL       -1

using boolean = bool;

class pnt {
public:
   double x   = 0.0;              /* x coordinate */
   double y   = 0.0;              /* y coordinate */
   int id     = NIL;              /* original id of pnt */
   boolean in = false;            /* is in interior of CHs */
   friend std::ostream& operator<< (std::ostream& os, const pnt& p);
};

using Pnts = std::vector<pnt>;

typedef struct {
   int nde;
   int num;
} loop;


typedef struct {
   int vtx;               /* index in vtx[]/pnts[]  */
   int next;              /* next CCW node          */
   int prev;              /* next CW node           */ 
} node;


class rt_options {
public:
   boolean verbose         	= false;
   boolean index		   	= false;
   long seed				= NIL;
   int counter				= 1;
   int timeout				= 0;
   std::string input_file   = "";
   std::string output_file  = "";
   boolean randomized		= true;
   boolean onion			= false;
   boolean obj				= false;

   boolean timings			= false;

   boolean use_stdin		= false;

   int partition			= 1;
   boolean flip_tris		= false;
};

class pCmpX {
public:
  bool operator() (pnt a, pnt b) {
	  if      (a.x < b.x)       return true;
	    else if (a.x > b.x)     return false;
	    else  {
	       if      (a.y < b.y)  return true;
	       else if (a.y > b.y)  return false;
	       else                 return false;
	    }
  }
};
class pCmpY {
public:
  bool operator() (pnt a, pnt b) {
	  if      (a.y < b.y)       return true;
	    else if (a.y > b.y)     return false;
	    else  {
	       if      (a.x < b.x)  return true;
	       else if (a.x > b.x)  return false;
	       else                 return false;
	    }
  }
};

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

#endif
