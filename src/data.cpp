#include "data.h"

#include <set>
#include <algorithm>    // std::for_each

std::ostream& operator<< (std::ostream& os, const ExtremPnts& ep) {
	os << ep.xMinIdx << ", " << ep.xMaxIdx << ", " << ep.yMinIdx << ", " << ep.yMaxIdx;

	return os;
}


void Broker::partition(int num_sets) {
	assert(num_sets > 1);

	int pnts_per_set = num_pnts/num_sets + 1;

	/* out of pure laziness we start by only partition the x-values */
	for(int i = 0; i < num_sets; ++i) {
		auto set = Data(pnts_per_set);
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

	for(unsigned long i = 0; i+1 < sets.size(); ++i) {
		if(cfg->verbose) {
			std::cout << "merging " << i << " and " << i+1 << std::endl;
		}

		mergeSets(i,i+1);
		repairOnionZeroAfterMerge(i,i+1);
	}
}

void Broker::mergeSets(const uint i, const uint j) {
	auto setA = &sets[i]; auto setB = &sets[j];

	auto itA = setA->onionZero.begin();
	auto itB = setB->onionZero.begin();
	auto itA2 = itA, itAend = itA;
	auto itB2 = itB, itBend = itB;

	std::advance(itA,setA->ep.yMaxIdx);
	std::advance(itB,setB->ep.yMaxIdx);

	std::advance(itAend,setA->ep.yMinIdx);
	std::advance(itBend,setB->ep.yMinIdx);

	do {
		itB2 = std::next(itB);
		if(itB2 == setB->onionZero.end()) {
			itB2 = setB->onionZero.begin();
		}

		if(itA == setA->onionZero.begin()) {
			itA2 = std::prev(setA->onionZero.end());
		} else {
			itA2 = std::prev(itA);
		}

		pnt* pA  = &setA->pnts[itA->vtx];
		pnt* pB  = &setB->pnts[itB->vtx];
		pnt* pA2 = &setA->pnts[itA2->vtx];
		pnt* pB2 = &setB->pnts[itB2->vtx];

		bool cA  = CCW(pB,pA,pA2);	bool cB  = CCW(pB2,pB,pA);
		bool cA2 = CCW(pA,pA2,pB2); bool cB2 = CCW(pA2,pB2,pB);

		if(fourConvexPoints(pA,pA2,pB2,pB)) {
			/* both CCW, add quad face :) */
			Face f = {pA->id,pA2->id,pB2->id,pB->id};
			faces.push_back(f);
			/* iterate both */
			itA = itA2;
			itB = itB2;
		} else if(cA && (cB2 || pA2->y > pB2->y)) {
			/* chose B,A,A2 */
			Face f1 = {pB->id,pA->id,pA2->id};
			faces.push_back(f1);
			/* iterate single */
			itA = itA2;
		} else if( cB && (cA2 || pA2->y < pB2->y ) ) {
			/* chose B2,B,A */
			Face f1 = {pB2->id,pB->id,pA->id};
			faces.push_back(f1);
			/* iterate single */
			itB = itB2;
		} else {
			/* one (initial) CH point may be far above/below other */
			if(pA->y > pB->y) {
				itA = itA2;
			} else {
				itB = itB2;
			}
			if(cfg->verbose) {
				std::cout << "ERROR! some bad decisions were made... (mergeSets)" << std::endl;
			}
		}
	} while(itA != itAend && itB != itBend);

	/* add the final triangle */
	std::set<pnt*> triPnts = {&setA->pnts[itA->vtx],
			&setB->pnts[itB->vtx],
			&setA->pnts[itA2->vtx],
			&setB->pnts[itB2->vtx]
	};

	Face f;
	for(auto p : triPnts) {f.push_back(p->id);}
	faces.push_back(f);
}

void Broker::writeFacesToFile(FILE *output) const {
	std::stringstream ss;
	for(auto f : faces) {
		ss << "f";
		for(auto i : f) {
			ss << " " << i+1;
		}
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

void Broker::repairOnionZeroAfterMerge(const uint i, const uint j) {
	auto setA = &sets[i]; auto setB = &sets[j];

	auto itA = setA->onionZero.begin();
	auto itB = setB->onionZero.begin();
	auto itAend = itA;
	auto itBend = itB;

	std::advance(itA,setA->ep.yMaxIdx);
	std::advance(itB,setB->ep.yMaxIdx);
	std::advance(itAend,setA->ep.yMinIdx);
	std::advance(itBend,setB->ep.yMinIdx);

	auto itAmerge = itA;
	auto itBmerge = itBend;

	std::cout << "begin removal..." << std::endl; fflush(stdout);

	/*			removing the merged part of the zero onion layer 		*/
	itA = setA->cyclicNext(itA);
	itB = setB->cyclicPrev(itB);
	itAend = setA->cyclicPrev(itAend);
	itBend = setB->cyclicNext(itBend);

	if(itAend - itA > 0) {
		std::cout << "a"; fflush(stdout);
		setA->onionZero.erase(itA,itAend);
	} else {
		std::cout << "a1"; fflush(stdout);
		setA->onionZero.erase(itAend,setA->onionZero.end());
		setA->onionZero.erase(setA->onionZero.begin(),itA);
	}
	if(itBend - itB > 0) {
		std::cout << "b"; fflush(stdout);
		setB->onionZero.erase(itB,itBend);
	} else {
		std::cout << "b1"; fflush(stdout);
		setB->onionZero.erase(itBend,setB->onionZero.end());
		setB->onionZero.erase(setB->onionZero.begin(),itB);
	}

	std::cout << "begin replacement..." << std::endl; fflush(stdout);

	/*			repair the zero onion layer of both sets		 		*/
	itA = setA->onionZero.begin();
	itB = setB->onionZero.begin();
	itAend = itA; itBend = itB;
	std::advance(itA,setA->ep.yMaxIdx);
	std::advance(itB,setB->ep.yMaxIdx);
	std::advance(itAend,setA->ep.yMinIdx);
	std::advance(itBend,setB->ep.yMinIdx);

	/* backup one list */
//	std::list<node> listBak(setA->onionZero);

	std::cout << "A..." << std::endl; fflush(stdout);
	if(itBend - itB > 0) {
		setA->onionZero.insert(itAend,itB,itBend);
	} else {
		setA->onionZero.insert(itAend,itB,setB->onionZero.end());
		setA->onionZero.insert(itAend,setB->onionZero.begin(),itBend);
	}
	std::cout << "B..." << std::endl; fflush(stdout);

	if(itAend - itA > 0) {
		setB->onionZero.insert(itBend,itA,itAend);
	} else {
		setB->onionZero.insert(itBend,itA,setA->onionZero.end());
		setB->onionZero.insert(itBend,setA->onionZero.begin(),itAend);
	}
	std::cout << "end rep. ..." << std::endl; fflush(stdout);
}

void Data::backupOnionZero(int idx) {
	auto l = layers[idx];
	node startNode = nodes[l.nde];
	node nIt = startNode;
	do {
		onionZero.push_back(nIt);
		nIt = nodes[nIt.next];
	} while(nIt.vtx != startNode.vtx);
}


void Data::determineExtremPoints() {
	auto minX = pnts[0].x;	auto maxX = pnts[0].x;
	auto minY = pnts[0].y;	auto maxY = pnts[0].y;

	int idx = 0;
	for(auto n : onionZero) {
		auto p = pnts[n.vtx];
		if(p.x < minX) {minX = p.x; ep.xMinIdx = idx;}
		if(p.x > maxX) {maxX = p.x; ep.xMaxIdx = idx;}
		if(p.y < minY) {minY = p.y; ep.yMinIdx = idx;}
		if(p.y > maxY) {maxY = p.y; ep.yMaxIdx = idx;}
		++idx;
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
		auto p = pnts[n.vtx];
		std::cout << " " <<  p.id+1;
	}
	std::cout << std::endl;
}

