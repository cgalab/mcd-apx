#include "data.h"

#include <set>
#include <algorithm>    // std::for_each

std::ostream& operator<< (std::ostream& os, const ExtremPnts& ep) {
//	os << ep.xMinIdx << ", " << ep.xMaxIdx << ", " << ep.yMinIdx << ", " << ep.yMaxIdx;

	return os;
}


void Broker::partition(int num_sets) {
	assert(num_sets > 1);

	int pnts_per_set = num_pnts/num_sets + 1;

	/* out of pure laziness we start by only partition the x-values */
	for(int i = 0; i < num_sets; ++i) {
		auto set = Data(pnts_per_set);
		if(i == 0) {set = Data();}
		for(int j = 0; j < pnts_per_set && (i*pnts_per_set + j) < num_pnts; ++j) {
			set.pnts[j] = pnts[i*pnts_per_set + j];
			set.num_pnts = j;
		}
		++set.num_pnts;
		sets.push_back(set);
	}
}

void Broker::merge() {
	for(auto& s : sets) {
		s.determineExtremPoints();
		if(cfg->verbose) {
			s.printLayer();
		}
	}

	if(cfg->verbose) {
		for(auto& s : sets) {
			std::cout << s.ep << std::endl;
		}
	}

	for(unsigned long i = 1; i < sets.size(); ++i) {
		if(cfg->verbose) {
			std::cout << "merging " << 0 << " and " << i << std::endl;
		}

//		mergeSets(0,i);
	}
}

void Broker::mergeSets(const uint i, const uint j) {

}

void Broker::writeFacesToFile(FILE *output) const {
	std::stringstream ss;
	for(auto f : faces) {
		ss << "f";
		for(auto i : f) {ss << " " << i+1;}
		ss << "\n";
	}
	fprintf(output, "%s", ss.str().c_str());
}


void Broker::printSets() const {
	for(auto s : sets) {
		s.printPnts();
	}
}


bool Broker::fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *pd) {
	bool cw  = ( CW(pa,pb,pc) &&  CW(pb,pc,pd) &&  CW(pc,pd,pa) &&  CW(pd,pa,pb));
	bool ccw = (CCW(pa,pb,pc) && CCW(pb,pc,pd) && CCW(pc,pd,pa) && CCW(pd,pa,pb));
	return cw || ccw;
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


void Data::backupOnionZero(int idx) {
	auto l = layers[idx];
	node startNode = nodes[l.nde];
	node nIt = startNode;
	do {
		onionZero.push_back(nIt.vtx);
		nIt = nodes[nIt.next];
	} while(nIt.vtx != startNode.vtx);
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

