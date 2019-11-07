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

	mergeSomeTris();

	for(auto pair : mergeTriToFace) {
		auto ta = tri.triangles[pair[0]];
		auto tb = tri.triangles[pair[1]];
		auto comPair = tri.getCommonPair(ta,tb);
		int aIdx = tri.getMissingCorner(ta,comPair[0],comPair[1]);
		int cIdx = tri.getMissingCorner(tb,comPair[0],comPair[1]);

		faces.push_back(Face({ aIdx, comPair[0], cIdx, comPair[1] }));
	}

	for(unsigned long i = 0; i < tri.triangles.size(); ++i) {
		if(visitedTris.find(i) == visitedTris.end()) {
			auto t = tri.triangles[i];
			faces.push_back(Face({t.a,t.b,t.c}));
		}

	}
}


void Broker::mergeSomeTris() {
	std::vector<unsigned long> triQueue;
	for(unsigned long i=0; i < tri.triangles.size();++i) {triQueue.push_back(i);}


	//std::random_shuffle ( triQueue.begin(), triQueue.end(), myrandom);
	auto rng = std::default_random_engine {};
	std::shuffle(std::begin(triQueue), std::end(triQueue), rng);

	for(auto idx : triQueue) {
		if(visitedTris.find(idx) == visitedTris.end()) {
			attemptMerge(idx);
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

void Broker::attemptMerge(int triIdx) {
	auto tringles = tri.triangles;
	const Triangle t = tringles[triIdx];
	std::vector<int> n;
	for(auto idx : {t.nAB,t.nBC,t.nCA}) {
		if(idx != NIL && visitedTris.find(idx) == visitedTris.end()) {
			n.push_back(idx);
		}
	}
	for(auto idx : n) {
		const Triangle tn = tringles[idx];
		if(tri.isConvexQuad(t,tn))  {
			visitedTris.insert(triIdx);
			visitedTris.insert(idx);
			mergeTriToFace.push_back({triIdx,idx});
			return;
		}
	}
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
