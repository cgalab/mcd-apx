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

using OnionZero    = std::list<int>;
using NodeIterator = OnionZero::iterator;
using Onions = std::vector<OnionZero>;

struct ExtremPnts {
	/* onionZero Iterator */
	NodeIterator xMinIdx, xMaxIdx, yMinIdx, yMaxIdx;
	friend std::ostream& operator<< (std::ostream& os, const ExtremPnts& ep);
};

int p_comp(const void *, const void *);

class Data {
public:
	Data(const Pnts &pntsVector);
	Data(const Data &data);
	Data(int size = 0);
	~Data();

	Data& operator=(const Data& other);

	pnt getInnerPnt();

	int num_pnts    = 0;
	int num_layers  = 0;
	int num_nodes   = 0;
	int lower_bound = 0;

	void printPnts() const;
	void printLayer() const;

	void backupOnionZero(int idx = 0);
	void determineExtremPoints();

	void resortPntsX() {qsort(pnts, num_pnts, sizeof(pnt), p_comp);};
	OnionZero getZeroWithOritinalID();

	pnt  *pnts;
	loop *layers;
	node *nodes;

	ExtremPnts ep;

	OnionZero onionZero;

	NodeIterator cyclicNext(NodeIterator it);
	NodeIterator cyclicPrev(NodeIterator it);
};


#endif
