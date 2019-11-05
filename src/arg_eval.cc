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

#include <string>       // std::string
#include <iostream>     // std::cout

/*                                                                           */
/* get standard libraries                                                    */
/*                                                                           */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

/*                                                                           */
/* get my header files                                                       */
/*                                                                           */
#include "defs.h"
#include "headers.h"


using Args = std::pair<int,char**>;

static struct option long_options[] = {
		{ "help"        	, no_argument      , 0, 'h'},
		{ "output"          , required_argument, 0, 'o'},
		{ "input"   	    , required_argument, 0, 'i'},
		{ "obj"   	    	, no_argument      , 0, 'w'},
		{ "onion"   	    , no_argument      , 0, 'O'},
		{ "seed"   	    	, required_argument, 0, 's'},
		{ "counter"   	    , required_argument, 0, 'c'},
		{ "timeout"   	    , required_argument, 0, 'T'},
		{ "index"   	    , no_argument      , 0, 'I'},
		{ "random"   	    , no_argument      , 0, 'R'},
		{ "verbose" 	    , no_argument      , 0, 'v'},
		{ "timings"     	, no_argument      , 0, 't'},
		{ 0, 0, 0, 0}
};

[[noreturn]]
 static void
 usage(const char *progname, int err) {
	FILE *f = err ? stderr : stdout;

	fprintf(f,"Usage: %s [options] <input file>\n", progname);
	fprintf(f,"  Options: --output  \t| --o <filename> \t write output\n");
	fprintf(f,"           --input   \t| --i <file>\t\t (option) since the input can also be simple added as last parameter\n");
	fprintf(f,"           --obj     \t| --w \t\t\t write output in wavefront obj format\n");
	fprintf(f,"           --onion   \t| --O \t\t\t use onion approach\n");
	fprintf(f,"           --seed    \t| --s NUM\t\t define fixed seed for random approach\n");
	fprintf(f,"           --counter \t| --c NUM\t\t set a counter (default: 1)\n");
	fprintf(f,"           --timeout \t| --T NUM\t\t set a timeout\n");
	fprintf(f,"           --index   \t| --I \t\t\t enable index\n");
	fprintf(f,"           --random  \t| --R \t\t\t use randomized approach (default)\n");
	fprintf(f,"           --verbose \t| --v \t\t\t print processing information\n");
	fprintf(f,"           --timings \t| --t \t\t\t print timings [ms]\n");
	fprintf(f,"\n");
	fprintf(f,"Input format is a list of coordinates, one pair of x y in each line.\n");
	fprintf(f,"\n");
	exit(err);
}

/*                                                                           */
/* this function parses the command-line arguments                           */
/*                                                                           */
void ArgEval(int argc, char *argv[], rt_options *rt_opt)
{

//   rt_opt->seed    = NIL;
//   rt_opt->counter = 1;
//   rt_opt->timeout = 0;
//   rt_opt->verbose = false;
//   rt_opt->index   = false;
//
//   rt_opt->randomized = true;
//   rt_opt->onion = false;

	while (1) {
		int option_index = 0;
		int r = getopt_long(argc, argv, "hO:R:", long_options, &option_index);

		if (r == -1) break;
		switch (r) {
		case 'h':
			usage(argv[0], 0);
			break;

		case 'o':
			rt_opt->output_file = std::string(optarg);
			break;

		case 'i':
			rt_opt->input_file = std::string(optarg);
			break;

		case 'w':
			rt_opt->obj = true;
			break;

		case 'O':
			rt_opt->onion 		= true;
			rt_opt->randomized  = false;
			break;

		case 's':
			rt_opt->seed = atoi(optarg);
			break;

		case 'c':
			rt_opt->counter = atoi(optarg);
			break;

		case 'T':
			rt_opt->timeout = atoi(optarg);
			break;

		case 'I':
			rt_opt->index = true;
			break;

		case 'R':
			rt_opt->onion 		= false;
			rt_opt->randomized  = true;
			break;

		case 'v':
			rt_opt->verbose = true;
			break;


		case 't':
			rt_opt->timings = true;
			std::cout << "not implemented." << std::endl;
			break;

		default:
			std::cerr << "Invalid option " << (char)r << std::endl;
			exit(1);
		}
	}

	if(rt_opt->input_file != "") {return;}

	if (argc - optind > 1) {usage(argv[0], 1);}

	if (argc - optind == 1) {
		std::string fn(argv[optind]);
		if (fn != "-") {
			rt_opt->input_file = fn;
		}
	}

   return;
}
