#include "data.h"

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

	for(int i = 0; i+1 < sets.size(); i+=2) {
		if(cfg->verbose) {
			std::cout << "merging " << i << " and " << i+1 << std::endl;
		}

		mergeSets(i,i+1);
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

		if(fourConvexPoints(pA,pA2,pB2,pB)) {
			/* both CCW, add quad face :) */
			Face f = {pA->id,pA2->id,pB2->id,pB->id};
			faces.push_back(f);
			/* iterate both */
			if(itA2 != itAend) {itA = itA2;}
			if(itB2 != itBend) {itB = itB2;}
		} else {
			/* at least one corner is 'reflex' */
			bool cA  = CCW(pB,pA,pA2);	bool cB  = CCW(pB2,pB,pA);
			bool cA2 = CCW(pA,pA2,pB);	bool cB2 = CCW(pA2,pB2,pB);

			if(cA || cB) {
				if((cA && cB && pA2->y > pB2->y) || (cA2 && !cB2)) {
					/* chose B,A,A2 */
					Face f1 = {pB->id,pA->id,pA2->id};
					faces.push_back(f1);
					/* iterate single */
					if(itA2 != itAend) {itA = itA2;}
				} else if( (cA2 && cB2 && pA2->y < pB2->y) || (!cA2 && cB2)) {
					/* chose B2,B,A */
					Face f1 = {pB2->id,pB->id,pA->id};
					faces.push_back(f1);
					/* iterate single */
					if(itB2 != itBend) {itB = itB2;}
				}
			}

			if(CCW(pA,pA2,pB2) && CCW(pB2,pB,pA) ) {
				Face f1 = {pA->id,pA2->id,pB2->id};
				Face f2 = {pB2->id,pB->id,pA->id};
				faces.push_back(f1);
				faces.push_back(f2);
				/* iterate both */
				if(itA2 != itAend) {itA = itA2;}
				if(itB2 != itBend) {itB = itB2;}
			} else if(CCW(pA2,pB2,pB) && CCW(pB,pA,pA2)) {
				Face f1 = {pA2->id,pB2->id,pB->id};
				Face f2 = {pB->id,pA->id,pA2->id};
				faces.push_back(f1);
				faces.push_back(f2);
				/* iterate both */
				if(itA2 != itAend) {itA = itA2;}
				if(itB2 != itBend) {itB = itB2;}
			} else {
				if(cfg->verbose) {
					std::cout << "ERROR! some bad decisions were made..." << std::endl;
//					throw UNKNOWN_ERROR;
				}
			}
		}
	} while(itA2 != itAend && itB2 != itBend);
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

void Data::backupOnionZoro(int idx) {
	auto l = layers[idx];
	node startNode = nodes[l.nde];
	node nIt = startNode;
	do {
		onionZero.push_back(nIt);
		nIt = nodes[nIt.next];
	} while(nIt.vtx != startNode.vtx);
}

void Data::printPnts() const {
	std::cout << "num_pnts = "    << num_pnts    << std::endl
			  << "num_layers = "  << num_layers  << std::endl
			  << "num_nodes   = " << num_nodes   << std::endl
			  << "lower_bound = " << lower_bound << std::endl;
	for(int i = 0; i < num_pnts; ++i) {
		std::cout
			<< "id: " << pnts[i].id
			<< " x: " << pnts[i].x
			<< " y: " << pnts[i].y << std::endl;
	}
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


void Data::printLayer() const {
	std::cout << "f";
	for(auto n : onionZero) {
		auto p = pnts[n.vtx];
		std::cout << " " <<  p.id+1;
	}
	std::cout << std::endl;
}



