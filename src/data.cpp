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

		mergeSets(0,i);
	}
}

void Broker::mergeSets(const uint i, const uint j) {
	auto setA = &sets[i]; auto setB = &sets[j];

	auto itA = setA->ep.yMaxIdx;
	auto itB = setB->ep.yMaxIdx;

	auto itA2 = itA;
	auto itB2 = itB;

	auto itAend = setA->ep.yMinIdx;
	auto itBend = setB->ep.yMinIdx;

	std::vector<NodeIterator> RemoveFromA;
	std::vector<NodeIterator> RemoveFromB;

	bool A2done = false, B2done = false;

	do {
		if(itA2 != itAend) {itA2 = setA->cyclicPrev(itA);}
		if(itB2 != itBend) {itB2 = setB->cyclicNext(itB);}

		pnt* pA  = &setA->pnts[itA->vtx];
		pnt* pB  = &setB->pnts[itB->vtx];
		pnt* pA2 = &setA->pnts[itA2->vtx];
		pnt* pB2 = &setB->pnts[itB2->vtx];

		bool cA  = CCW(pB,pA,pA2);	bool cB  = CCW(pB2,pB,pA);
		bool cA2 = CCW(pA,pA2,pB2); bool cB2 = CCW(pA2,pB2,pB);
		if((!A2done || itA2 != itAend) && (!B2done || itB2 != itBend) && fourConvexPoints(pA,pA2,pB2,pB)) {
			/* both CCW, add quad face :) */
			Face f = {pA->id,pA2->id,pB2->id,pB->id};
			faces.push_back(f);
			/* iterate both */
			if(itA2 == itAend) {A2done = true;}
			if(itB2 == itBend) {B2done = true;}
			if(itA2 != itAend) {itA = itA2; RemoveFromA.push_back(itA2);}
			if(itB2 != itBend) {itB = itB2; RemoveFromB.push_back(itB2);}
		} else if( (!A2done || itA2 != itAend) && cA && !PntInTri(pB2,pA2,pA,pB)) { //(cB2 || pA2->y > pB2->y || pB2->x > pB->x)) {
			/* chose B,A,A2 */
			Face f1 = {pB->id,pA->id,pA2->id};
			faces.push_back(f1);
			/* iterate single */
			if(itA2 == itAend) {A2done = true;}
			if(itA2 != itAend) {itA = itA2; RemoveFromA.push_back(itA2);}
		} else if((!B2done || itB2 != itBend) && cB && !PntInTri(pA2,pA,pB,pB2)) { //(cA2 || pA2->y < pB2->y || pA2->x < pA->x ) ) {
			/* chose B2,B,A */
			Face f1 = {pB2->id,pB->id,pA->id};
			faces.push_back(f1);
			/* iterate single */
			if(itB2 == itBend) {B2done = true;}
			if(itB2 != itBend) {itB = itB2; RemoveFromB.push_back(itB2);}
		} else {
			/* one (initial) CH point may be far above/below other */
			if(itA2 != itAend && CW(pA,pA2,pB)) {
				std::cout<< " it A ";
				itA = itA2;
			} else if(itB2 != itBend && CCW(pB,pB2,pA)) {
				std::cout<< " it B ";
				itB = itB2;
			} else if(cfg->verbose && itA2 != itAend && itB2 != itBend) {
				std::cout << "ERROR! (" << i << "," << j << ") some bad decisions were made... (mergeSets)" << std::endl;
			}
		}
	} while(!(itA2 == itAend && itB2 == itBend));

	std::cout << std::boolalpha << "a2: " << A2done << ", b2:" << B2done << std::endl;

	/* add the final triangle */
	std::vector<pnt*> pF = {
			&setA->pnts[itA2->vtx],
			&setB->pnts[itB2->vtx]
	};

	if(!B2done) {
		pF.push_back(&setB->pnts[itB->vtx]);
	}
	if(!A2done) {
		pF.push_back(&setA->pnts[itA->vtx]);
	}

	std::set<pnt*> triSet;
	std::cout << " points: "<<std::endl;
	std::for_each(pF.begin(),pF.end(),[&](pnt *p) {std::cout<< *p << std::endl; triSet.insert(p);});

	if(triSet.size() == 4 && fourConvexPoints(pF[0],pF[1],pF[2],pF[3])) {
		std::cout << "-add a convex 4-face! " << std::endl;
		faces.push_back({pF[0]->id,pF[1]->id,pF[2]->id,pF[3]->id});
	} else if(triSet.size() == 4) {
		int cwIdx = NIL;
		for(int i = 0; i < 4; ++i) {
			if(CW(pF[i],pF[(i+1)%4],pF[(i+2)%4])) {
				cwIdx = (i+1)%4;
				break;
			}
		}
		if(cwIdx == NIL) {std::cout << "ERROR cw idx not found!" << std::endl;}

		faces.push_back({pF[(cwIdx)]->id,pF[(cwIdx+1)%4]->id,pF[(cwIdx+2)%4]->id});
		faces.push_back({pF[(cwIdx+2)%4]->id,pF[(cwIdx+3)%4]->id,pF[(cwIdx+4)%4]->id});
		std::cout << " add two faces ! " << std::endl;
	} else {
		Face f;
		std::cout << std::endl << " add a single face ! ";
		std::for_each(triSet.begin(),triSet.end(),[&](pnt *p) {std::cout << p->id << " "; f.push_back(p->id);});
//		faces.push_back({pF[0]->id,pF[1]->id,pF[2]->id});
		faces.push_back(f);
	}


	/* ----------------------------- update onionZero --------------------------- */

	/* remove nodes handled in the merge above */
	std::for_each(RemoveFromA.begin(),RemoveFromA.end(),[&](NodeIterator i){setA->onionZero.erase(i);});
	std::for_each(RemoveFromB.begin(),RemoveFromB.end(),[&](NodeIterator i){setB->onionZero.erase(i);});

	/* itAend is the start of our insert of the list that starts at itBend */
	itB = itBend; itA = itAend;
	auto itBStop = itB;
	do {
		itA = setA->cyclicNext(itA);
		/* add point */
		auto pnt = setB->pnts[itB->vtx];
		setA->pnts[setA->num_pnts] = pnt;
		node n = {setA->num_pnts,NIL,NIL};
		++setA->num_pnts;
		itA = setA->onionZero.insert(itA,n);

		itB = setB->cyclicNext(itB);
	} while(itB != itBStop);

	std::cout << "A onion: "; setA->printLayer(); std::cout << std::endl;

	setA->determineExtremPoints();
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

	for(auto it = onionZero.begin(); it != onionZero.end(); ++it) {
		auto p = pnts[it->vtx];
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
		auto p = pnts[n.vtx];
		std::cout << " " <<  p.id+1;
	}
	std::cout << std::endl;
}

