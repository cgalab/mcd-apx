#ifndef BROKER_H_
#define BROKER_H_

#include <assert.h>
#include <iostream>
#include <cstdlib>

#include <list>
#include <vector>
#include <sstream>
#include <string>

#include "defs.h"
#include "numerics.h"
#include "data.h"
#include "Tri.h"


class Broker : public Data {
	using Sets  = std::vector<Data>;
	using Cfg   = rt_options;
	using Face  = std::vector<long>;
	using Faces = std::vector<Face>;

public:
	Broker(int size = 0):Data(size) {}

	void partition(int num_sets = 1);
	void collectZeroOnions();
	Pnts collectHolePnts();

	void merge();
	void writeFacesToFile(FILE *output) const;

	void printSets() const;

	void mergeSomeTris();

	Sets  sets;
	Faces faces; /* faces from the merge */

	Cfg*  cfg = nullptr;

	bool fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *d);

	Tri tri;

private:
	std::set<int> visitedTris;
	std::vector< std::vector<int>>  mergeTriToFace;

	Onions allZeroOnions;

	void attemptMerge(int triIdx);

//	int myrandom (int i) { return std::rand()%i;}
};

#endif
