#ifndef DATA_H_
#define DATA_H_

#include <assert.h>
#include <iostream>

#include <list>
#include <vector>
#include <sstream>
#include <string>

#include "defs.h"
#include "numerics.h"

struct ExtremPnts {
	int xMinIdx, xMaxIdx, yMinIdx, yMaxIdx;
	friend std::ostream& operator<< (std::ostream& os, const ExtremPnts& ep);
};

class Data {
public:
	Data(const Data &data):Data(data.num_pnts) {
		num_pnts    = data.num_pnts;
		num_layers  = data.num_layers;
		num_nodes   = data.num_nodes;
		lower_bound = data.lower_bound;
		std::copy(&data.pnts[0],   &data.pnts[0]   + num_pnts, &pnts[0]);
		std::copy(&data.layers[0], &data.layers[0] + num_pnts, &layers[0]);
		std::copy(&data.nodes[0],  &data.nodes[0]  + num_pnts, &nodes[0]);
		onionZero = data.onionZero;
		ep = data.ep;
	}
	Data(int size = 0) {
		if(size == 0) {
			pnts   = new pnt[MAX];
			layers = new loop[MAX];
			nodes  = new node[MAX];
		} else {
			pnts   = new pnt[size];
			layers = new loop[size];
			nodes  = new node[size];
		}
	}
	~Data() {
		delete[] pnts;
		delete[] nodes;
		delete[] layers;
	}

	Data& operator=(const Data& other) {
		num_pnts = other.num_pnts;
		std::copy(&other.pnts[0],   &other.pnts[0]   + num_pnts, &pnts[0]);
		std::copy(&other.layers[0], &other.layers[0] + num_pnts, &layers[0]);
		std::copy(&other.nodes[0],  &other.nodes[0]  + num_pnts, &nodes[0]);
		return *this;
	}

	int num_pnts    = 0;
	int num_layers  = 0;
	int num_nodes   = 0;
	int lower_bound = 0;

	void printPnts() const;
	void printLayer() const;

	void backupOnionZoro(int idx = 0);
	void determineExtremPoints();

	pnt  *pnts;
	loop *layers;
	node *nodes;

	ExtremPnts ep = {0,0,0,0};

	std::vector<node> onionZero;
};


class Broker : public Data {
	using Sets  = std::vector<Data>;
	using Cfg   = rt_options;
	using Face  = std::vector<int>;
	using Faces = std::vector<Face>;

public:
	Broker(int size = 0):Data(size) {}

	void partition(int num_sets = 1);
	void merge();
	void writeFacesToFile(FILE *output) const;

	void printSets() const;

	Sets  sets;
	Faces faces; /* faces from the merge */

	Cfg*  cfg = nullptr;

	bool fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *d);

private:
	void mergeSets(const uint i, const uint j);
};

#endif
