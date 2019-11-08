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

	using FaceIterator = Face::iterator;

public:
	Broker(int size = 0):Data(size) {}

	void partition(int num_sets = 1);
	void merge();

	void mergeSomeTris();

	void runTriangleOnlyApproach();

	void collectZeroOnions();
	Pnts collectHolePnts();

	void writeFacesToFile(FILE *output) const;
	void printSets() const;


	Sets  sets;
	Faces faces; /* faces from the merge */

	Cfg*  cfg = nullptr;

	bool fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *d);

	Tri tri;

private:
	std::set<int> visitedTris;
	std::vector< std::vector<int>>  mergeTriToFace;

	Onions allZeroOnions;

	bool addTriToFace(long int tidx, Face& f);

	void attemptExpansion(int triIdx);

	FaceIterator cNext(Face& f, FaceIterator it) {
		return (++it == f.end()) ? f.begin() : it;
	}
	FaceIterator cPrev(Face& f, FaceIterator it) {
		return (it == f.begin()) ? (f.end()-1) : --it;
	}
};

#endif
