#include "broker.h"

#include <set>
#include <algorithm>    // std::for_each, std::sort
#include <math.h>
#include <random>

void Broker::partition(int num_sets) {
	assert(num_sets > 1);
	pCmpY pCmpY;

	unsigned long width = ceil(sqrt(num_sets));
	if(width < 2) {width = 2;}

	if(cfg->verbose) {
		std::cout << "using " << num_sets << " partitions: " << width << "x" << width << " grid" << std::endl;
	}

	int num_per_stripe = 1 + ceil(num_pnts/width);

	std::vector<Pnts> stripes;
	Pnts ptmp;
	for(int i = 0; i < num_pnts; ++i) {
		ptmp.push_back(pnts[i]);
		if((i%(num_per_stripe))+1 == num_per_stripe) {
			stripes.push_back(ptmp);
			ptmp.clear();
		}
	}


	if(!ptmp.empty()) {
		if(stripes.size() < width) {
			stripes.push_back(ptmp);
		} else {
			auto& lastStripe = stripes.back();
			for(auto p : ptmp) {
				lastStripe.push_back(p);
			}
		}
	}

	for(auto& s : stripes) {
		std::sort(s.begin(),s.end(), pCmpY);
	}

	ptmp.clear();
	unsigned long cnt = 0, row = 0;
	for(auto s : stripes) {
		++row;
		for(auto p : s) {
			ptmp.push_back(p);
			++cnt;
			if( cnt == num_pnts/(width*width)) {
				sets.push_back(Data(ptmp));
				ptmp.clear();
				cnt = 0;
			}
		}
		if(!ptmp.empty()) {
			if(sets.size() < row*width) {
				sets.push_back(Data(ptmp));
			} else {
				auto& lastData = sets.back();
				for(int i = 0 ; i < lastData.num_pnts; ++i) {
					ptmp.push_back(lastData.pnts[i]);
				}
				sets.pop_back();
				sets.push_back(Data(ptmp));
			}
			ptmp.clear();
		}
	}

	for(auto &s : sets) {
		s.resortPntsX();
	}
}

void Broker::merge() {

	collectZeroOnions();

	auto holePnts = collectHolePnts();

	tri.runTriangle(pnts,num_pnts,holePnts,allZeroOnions);

	std::srand(cfg->seed);

	if(cfg->flip_tris) {
		flipTriangles();
	}

	mergeSomeTris();

	for(unsigned long i = 0; i < tri.triangles.size(); ++i) {
		if(visitedTris.find(i) == visitedTris.end()) {
			auto t = tri.triangles[i];
			faces.push_back(Face({t.a,t.b,t.c}));
		}

	}
}

void Broker::flipTriangles() {

}

void Broker::mergeSomeTris() {
	std::vector<unsigned long> triQueue;
	for(unsigned long i=0; i < tri.triangles.size();++i) {triQueue.push_back(i);}


	auto rng = std::default_random_engine {};
	rng.seed(cfg->seed);
	std::shuffle(std::begin(triQueue), std::end(triQueue), rng);

	for(auto idx : triQueue) {
		if(visitedTris.find(idx) == visitedTris.end()) {
			attemptExpansion(idx);
		}
	}
}

Pnts Broker::collectHolePnts() {
	Pnts holePnts;

	for(auto s : sets) {
		holePnts.push_back(s.getInnerPnt());
	}
	return holePnts;
}

void Broker::collectZeroOnions() {
	for(auto& s : sets) {
		allZeroOnions.push_back(s.getZeroWithOritinalID());
	}
}

void Broker::attemptExpansion(int triIdx) {
	auto tringles = tri.triangles;
	const Triangle t = tringles[triIdx];

	/* we try to expand this face */
	Face f = {t.a,t.b,t.c};

	std::set<long int> checked;
	std::list<long int> candidates = {t.nAB,t.nBC,t.nCA};

	checked.insert(triIdx);

	do {
		auto cndIdx = candidates.front();
		if(cndIdx != NIL
			&& visitedTris.find(cndIdx) == visitedTris.end()
			&&     checked.find(cndIdx) == checked.end()
		) {

			/* check if we can add this tri and f stays convex */
			if(addTriToFace(cndIdx,f)) {

				/* we can add it, then we should check its neighbours as well */
				auto tCnd = tringles[cndIdx];
				for(auto tCndN : {tCnd.nAB,tCnd.nBC,tCnd.nCA}) {
					if(tCndN != NIL
						&& visitedTris.find(tCndN) == visitedTris.end()
						&&     checked.find(tCndN) == checked.end())
					{
						candidates.push_back(tCndN);
					}
				}
				/* add to never touch it again */
				visitedTris.insert(cndIdx);
			}
			/* add to locally checked tris */
			checked.insert(cndIdx);
		}

		candidates.pop_front();
	} while(!candidates.empty());

	if(f.size() > 3) {
		visitedTris.insert(triIdx);
		faces.push_back(f);
	}
}

bool Broker::addTriToFace(long int tidx, Face& f) {
	auto t = tri.triangles[tidx];
	/* finding the indices */
	auto itB = f.begin();
	while(itB != f.end() && !t.hasIndex(*itB)) {++itB;}

	auto itX = cNext(f,itB);
	if(t.hasIndex(*itX)) {

	} else {
		itX = cPrev(f,itB);
		if(t.hasIndex(*itX)) {
			itB = itX;
		}
	}

	if(itB == f.end()) {
		return false;
	}

	auto itA = cPrev(f,itB);
	auto itD = cNext(f,itB);
	auto itE = cNext(f,itD);

	auto idxC = t.getThirdIndex(*itB,*itD);

	if(idxC == NIL) {
		return false;
	}

	auto Pa = pnts[*itA];
	auto Pb = pnts[*itB];
	auto Pc = pnts[idxC];
	auto Pd = pnts[*itD];
	auto Pe = pnts[*itE];

	if(CCW(&Pa,&Pb,&Pc) && CCW(&Pc,&Pd,&Pe)) {
		f.insert(itD,idxC);
		return true;
	}

	return false;
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
	int cnt = 0;
	for(auto s : sets) {
		std::cout << ++cnt << " set: " << std::endl;
		s.printPnts();
		std::cout << std::endl;
	}
}

bool Broker::fourConvexPoints(pnt *pa, pnt *pb, pnt *pc, pnt *pd) {
	bool cw  = ( CW(pa,pb,pc) &&  CW(pb,pc,pd) &&  CW(pc,pd,pa) &&  CW(pd,pa,pb));
	bool ccw = (CCW(pa,pb,pc) && CCW(pb,pc,pd) && CCW(pc,pd,pa) && CCW(pd,pa,pb));
	return cw || ccw;
}
