#ifndef BROKER_H_
#define BROKER_H_

#include <assert.h>
#include <iostream>
#include <cstdlib>

#include <list>
#include <vector>
#include <sstream>
#include <string>

#include <unordered_set>
#include <unordered_map>

#include <random>

#include "defs.h"
#include "numerics.h"
#include "data.h"
#include "Tri.h"

using Face  = std::list<long>;
using Faces = std::vector<Face>;
using FaceIterator = Face::iterator;

class Broker : public Data {
	using Sets  = std::vector<Data>;
	using Cfg   = rt_options;

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
	Tri tri;

	long getNumFaces() const {
		long num = 0;
		for(auto& f : faces) {
			if(f.size() > 1) {
				++num;
			}
		}
		return num;
	}

private:
	std::set<int> visitedTris;

	std::unordered_map<long int, long int> triToFaceMap;
	std::unordered_map<long int, std::list<long int>> faceToTriMap;

	std::default_random_engine rng = std::default_random_engine {};

	std::list<long> freeFaceSpace;

	Onions allZeroOnions;

	bool addTriToFace(long int tidx, Face& f);

	inline long addFace(Face f) {
		if(freeFaceSpace.empty()) {
			faces.push_back(f);
			return faces.size()-1;
		} else {
			auto idx = freeFaceSpace.front();
			freeFaceSpace.pop_front();
			faces[idx] = f;
			return idx;
		}
	}

	inline Face removeFace(long idx) {
		freeFaceSpace.push_back(idx);
		auto f = faces[idx];
		faces[idx] = {{}};
		return f;
	}



	void attemptExpansion(int triIdx, std::unordered_set<long int> allowedTris = {{}});

	void attemptFlipping(TriQueue &triQueue, unsigned long flips, bool inSet = false);

	std::unordered_set<long int> getTrisOfFaces(std::unordered_set<long int>& trifaces);
	std::unordered_set<long int> getFacesOfTriangles(TriQueue &tris);
	TriQueue selectNumTrisBFS(long int triIdx, unsigned long int num);

	inline FaceIterator cNext(Face& f, FaceIterator it) {return (++it == f.end()) ? f.begin() : it;}
	inline FaceIterator cPrev(Face& f, FaceIterator it) {return (it == f.begin()) ? (std::prev(f.end())) : --it;}

	void addTriAsFace(long int triIdx);
	bool fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *d);
};

#endif
