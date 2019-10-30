#ifndef DATA_H_
#define DATA_H_

#include <assert.h>
#include <iostream>

#include <list>

#include "defs.h"


class Data {
public:
	Data(const Data &data):Data(data.num_pnts) {
		num_pnts = data.num_pnts;
		std::copy(&data.pnts[0],   &data.pnts[0]   + num_pnts, &pnts[0]);
		std::copy(&data.layers[0], &data.layers[0] + num_pnts, &layers[0]);
		std::copy(&data.nodes[0],  &data.nodes[0]  + num_pnts, &nodes[0]);
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
	void printLayer(int idx) const;

	pnt  *pnts;
	loop *layers;
	node *nodes;
};


class Broker : public Data {
	using Sets  = std::list<Data>;
	using Cfg   = rt_options;

public:
	Broker(int size = 0):Data(size) {}

	void partition(int num_sets = 1);
	void merge();

	void printSets() const;

	Sets sets;
	Cfg*  cfg = nullptr;
};

#endif
