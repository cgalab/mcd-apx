#include "data.h"

#include <set>
#include <algorithm>    // std::for_each

std::ostream& operator<< (std::ostream& os, const ExtremPnts& ep) {
	os << *ep.xMinIdx << ", " << *ep.xMaxIdx << ", " << *ep.yMinIdx << ", " << *ep.yMaxIdx;

	return os;
}
Data& Data::operator=(const Data& other) {
		num_pnts = other.num_pnts;
		std::copy(&other.pnts[0],   &other.pnts[0]   + num_pnts, &pnts[0]);
		std::copy(&other.layers[0], &other.layers[0] + num_pnts, &layers[0]);
		std::copy(&other.nodes[0],  &other.nodes[0]  + num_pnts, &nodes[0]);
		return *this;
	}

NodeIterator Data::cyclicPrev(NodeIterator it) {
	if(it == onionZero.begin()) {
		return std::prev(onionZero.end());
	}
	return std::prev(it);
}


NodeIterator Data::cyclicNext(NodeIterator it) {
	auto nIt = std::next(it);
	return (nIt != onionZero.end()) ? nIt : onionZero.begin();
}

Data::Data(const Data &data):Data(data.num_pnts) {
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

Data::Data(const Pnts &pntsVector) {
	num_pnts = pntsVector.size();
	pnts   = new pnt[num_pnts]();
	layers = new loop[num_pnts]();
	nodes  = new node[num_pnts]();
	for(auto i=0; i < num_pnts; ++i) {
		pnts[i] = pntsVector[i];
	}
}

Data::Data(int size) {
	if(size == 0) {
		pnts   = new pnt[MAX]();
		layers = new loop[MAX]();
		nodes  = new node[MAX]();
	} else {
		pnts   = new pnt[size]();
		layers = new loop[size]();
		nodes  = new node[size]();
	}
}

Data::~Data() {
	delete[] pnts;
	delete[] nodes;
	delete[] layers;
}

void Data::backupOnionZero(int idx) {
	auto l = layers[idx];
	node startNode = nodes[l.nde];
	node nIt = startNode;
	do {
		onionZero.push_back(nIt.vtx);
		nIt = nodes[nIt.next];
	} while(nIt.vtx != startNode.vtx);
}

pnt Data::getInnerPnt() {
	auto it = onionZero.begin();
	auto Pa = pnts[*it];
	++it; ++it;
	auto Pb = pnts[*it];

	return midpoint(&Pa, &Pb);
}

void Data::determineExtremPoints() {
	auto minX = pnts[0].x;	auto maxX = pnts[0].x;
	auto minY = pnts[0].y;	auto maxY = pnts[0].y;

	for(auto it = onionZero.begin(); it != onionZero.end(); ++it) {
		auto p = pnts[*it];
		if(p.x < minX) {minX = p.x; ep.xMinIdx = it;}
		if(p.x > maxX) {maxX = p.x; ep.xMaxIdx = it;}
		if(p.y < minY) {minY = p.y; ep.yMinIdx = it;}
		if(p.y > maxY) {maxY = p.y; ep.yMaxIdx = it;}
	}
}

OnionZero Data::getZeroWithOritinalID() {
	OnionZero onion;
	for(auto n : onionZero) {
		onion.push_back(pnts[n].id);
	}
	return onion;
}

void Data::printPnts() const {
	std::cout << "num_pnts = "    << num_pnts    << std::endl
			  << "num_layers = "  << num_layers  << std::endl
			  << "num_nodes   = " << num_nodes   << std::endl
			  << "lower_bound = " << lower_bound << std::endl;
	for(int i = 0; i < num_pnts; ++i) {
		std::cout << pnts[i] << std::endl;
	}
}


void Data::printLayer() const {
	std::cout << "f";
	for(auto n : onionZero) {
		auto p = pnts[n];
		std::cout << " " <<  p.id+1;
	}
	std::cout << std::endl;
}

