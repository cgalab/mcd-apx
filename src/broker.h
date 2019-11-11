#ifndef BROKER_H_
#define BROKER_H_

#include <assert.h>
#include <iostream>
#include <cstdlib>

#include <list>
#include <vector>
#include <sstream>
#include <string>
#include <unordered_map>

#include <random>

#include "defs.h"
#include "numerics.h"
#include "data.h"
#include "Tri.h"


class Broker : public Data {
	using Sets  = std::vector<Data>;
	using Cfg   = rt_options;
	using Face  = std::vector<long>;
	using Faces = std::list<Face>;

	using FaceIterator = Face::iterator;

public:
	Broker(int size = 0):Data(size) {}

	void partition(int num_sets = 1);
	void merge();
	void startHoleRecursion();

	void mergeSomeTris();

	void runTriangleOnlyApproach() {merge();}

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

	std::unordered_map<long int, long int> triToFaceMap;
	std::unordered_map<long int, long int> faceToTriMap;

	std::default_random_engine rng = std::default_random_engine {};

	Onions allZeroOnions;

	bool addTriToFace(long int tidx, Face& f);


	void attemptExpansion(int triIdx);

	void attemptFlipping(TriQueue &triQueue);

	TriQueue selectNumTrisBFS(long int triIdx, long int num);

	inline FaceIterator cNext(Face& f, FaceIterator it) {return (++it == f.end()) ? f.begin() : it;}
	inline FaceIterator cPrev(Face& f, FaceIterator it) {return (it == f.begin()) ? (f.end()-1) : --it;}
};

#endif
