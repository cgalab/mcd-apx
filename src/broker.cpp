#include "broker.h"

#include <set>
#include <algorithm>    // std::for_each, std::sort
#include <math.h>


void Broker::addTriAsFace(long int triIdx) {
	auto t = tri.getTriangle(triIdx);
	auto faceIdx = addFace(Face({t.a,t.b,t.c}));

	auto it = triToFaceMap.find(triIdx);
	if(it != triToFaceMap.end()) {
		triToFaceMap.erase(it);
	}
	triToFaceMap.insert({{triIdx, faceIdx}});


	std::list<long int> onlyi;
	onlyi.push_back(triIdx);
	auto it2 = faceToTriMap.find(faceIdx);
	if(it2 != faceToTriMap.end()) {
		faceToTriMap.erase(it2);
	}
	faceToTriMap.insert({{faceIdx, onlyi }});
}

void Broker::merge() {
	/* only relevant if partition != 0 */
	collectZeroOnions();

	auto holePnts = collectHolePnts();

	tri.runTriangle(pnts,num_pnts,holePnts,allZeroOnions);

	mergeSomeTris();

	for(int i = 0; i < tri.getNumTriangles(); ++i) {
		if(visitedTris.find(i) == visitedTris.end()) {
			addTriAsFace(i);
		}
	}

	if(cfg->recurse_holes) {
		startHoleRecursion();
	}
}

void Broker::startHoleRecursion() {
	/* let us start by a BFS of a random triangle and select sqrt n tris */
	TriQueue triQueue(tri.getNumTriangles());
	std::iota(triQueue.begin(),triQueue.end(),0);
	std::shuffle(std::begin(triQueue), std::end(triQueue), rng);
	triQueue.resize((int)sqrt(triQueue.size()));


	std::cout << "TRI IDX: " << triQueue.back() << std::endl;

//	while(!triQueue.empty()) {
		Faces backup;
		long numFaces = getNumFaces();
		long numRemovedFaces = 0;

		long retries = 1;

		while(retries-- > 0) {
			auto selectedTris = selectNumTrisBFS(triQueue.back(),(long int)sqrt(tri.getNumTriangles()));

			/* let us find the faces that contain these triangles at the moment */
			auto facesOfTris = getFacesOfTriangles(selectedTris);

			/* let us obtain all triangles that form these faces */
			auto allTris = getTrisOfFaces(facesOfTris);

			/* let us backup and delete these faces */
			numRemovedFaces = facesOfTris.size();
			for(auto i : facesOfTris) {
				backup.push_back(removeFace(i));
			}

			/* new we try to find a better solution for this tri-subset */
			TriQueue allTrisVect;
			visitedTris.clear();
			for(auto t : allTris) {
				visitedTris.insert(t);
				allTrisVect.push_back(t);
			}
//			attemptFlipping(allTrisVect,1,true);

			visitedTris.clear();
			std::shuffle(std::begin(allTrisVect), std::end(allTrisVect), rng);

			for(auto idx : allTrisVect) {
				if(visitedTris.find(idx) == visitedTris.end()) {
					attemptExpansion(idx,allTris);
				}
			}

			for(auto idx : allTrisVect) {
				if(visitedTris.find(idx) == visitedTris.end()) {
					addTriAsFace(idx);
				}
			}

			long numFacesNew = getNumFaces();

			std::cout << "old: " << numFaces << ", removed: " << numRemovedFaces << ", new: " << numFacesNew << std::endl;

			if(numFacesNew > numFaces) {
				std::cout << "restore..." <<std::endl;
				long backupCnt = 0;
				/* we restore to the previous state */
				for(auto i : facesOfTris) {
					faces[i] = backup[backupCnt++];
				}
			} else {
				std::cout << "improve... " << numFacesNew <<std::endl;
				backup.clear();
				numFaces = numFacesNew;
			}
			faces.resize(numFaces);
			visitedTris.clear();
		}
//		triQueue.pop_back();
//	}
}

std::unordered_set<long int> Broker::getTrisOfFaces(std::unordered_set<long>& trifaces) {
	std::unordered_set<long int> tris;
	for(auto f : trifaces) {
		auto it = faceToTriMap.find(f);
		if(it != faceToTriMap.end()) {
			for(auto idx : it->second) {
				tris.insert(idx);
			}
		}
	}
	return tris;
}

std::unordered_set<long int> Broker::getFacesOfTriangles(TriQueue &tris) {
	std::unordered_set< long int > set;
	for(auto t : tris) {
		auto it = triToFaceMap.find(t);
		if(it != triToFaceMap.end()) {
			set.insert(it->second);
		}
	}
	return set;
}


TriQueue Broker::selectNumTrisBFS(long int triIdx, unsigned long int num) {
	TriQueue tq;
	auto t = tri.getTriangle(triIdx);
	std::set<long int> checked;

	tq.push_back(triIdx);
	checked.insert(triIdx);
	unsigned long it = 0;

	while(it < num) {
		t = tri.getTriangle(tq[it]);
		for(auto n : {t.nAB,t.nBC,t.nCA}) {
			if(n != NIL && checked.find(n) == checked.end()) {
				tq.push_back(n);
				checked.insert(n);
			}
		}
		++it;
	}

	return tq;
}

void Broker::mergeSomeTris() {
	TriQueue triQueue(tri.getNumTriangles());
	std::iota(triQueue.begin(),triQueue.end(),0);

	/* set a fixed seed if one is provided */
	if(cfg->seed != NIL) {
		rng.seed(cfg->seed);
	}


	if(cfg->flip_tris != 0) {
		auto flips = cfg->flip_tris;

		if(cfg->flip_tris == -1) {
			triQueue.resize((int)sqrt(tri.getNumTriangles()));
			flips = 1;
		}

		attemptFlipping(triQueue, flips);

		if(cfg->flip_tris == -1) {
			triQueue.clear();
			std::iota(triQueue.begin(),triQueue.end(),0);
		}
	}

	std::shuffle(std::begin(triQueue), std::end(triQueue), rng);

	for(auto idx : triQueue) {
		if(visitedTris.find(idx) == visitedTris.end()) {
			attemptExpansion(idx);
		}
	}
}

void Broker::attemptFlipping(TriQueue &triQueue, unsigned long flips, bool inSet) {
	while(flips-- > 0) {
		std::shuffle(std::begin(triQueue), std::end(triQueue), rng);
		for(auto idx : triQueue) {
			if(inSet || visitedTris.find(idx) == visitedTris.end()) {
				auto t = tri.getTriangle(idx);
				std::vector<long> nV = {t.nAB,t.nBC,t.nCA};
				std::shuffle(std::begin(nV), std::end(nV), rng);
				for(auto n : nV) {
					if(n != NIL &&
							((!inSet && visitedTris.find(n) == visitedTris.end())
							||
							(inSet && visitedTris.find(n) != visitedTris.end()))
					) {
						auto tN = tri.getTriangle(n);
						if(tri.isConvexQuad(t,tN)) {
							tri.flipPair(t,tN);
							if(!inSet) {
								visitedTris.insert(idx);
								visitedTris.insert(n);
							}
							break;
						}
					}
				}
			}
		}
		if(!inSet) {visitedTris.clear();}
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

void Broker::attemptExpansion(int triIdx, std::unordered_set<long int> allowedTris) {
	auto t = tri.getTriangle(triIdx);

	/* we try to expand this face */
	Face f = {t.a,t.b,t.c};

	std::set<long int>  checked;
	std::list<long int> candidates = {t.nAB,t.nBC,t.nCA};
	std::list<long int> trisInFace = {triIdx};

	checked.insert(triIdx);

	do {
		auto cndIdx = candidates.front();
		if(cndIdx != NIL
			&& visitedTris.find(cndIdx) == visitedTris.end()
			&&     checked.find(cndIdx) == checked.end()
			&& (allowedTris.empty() || allowedTris.find(cndIdx) != allowedTris.end())
		) {
			/* check if we can add this tri and f stays convex */
			if(addTriToFace(cndIdx,f)) {

				trisInFace.push_back(cndIdx);

				/* we can add it, then we should check its neighbours as well */
				auto tCnd = tri.getTriangle(cndIdx);
				for(auto tCndN : {tCnd.nAB,tCnd.nBC,tCnd.nCA}) {
					if(tCndN != NIL
						&& visitedTris.find(tCndN) == visitedTris.end()
						&&     checked.find(tCndN) == checked.end()
						&& (allowedTris.empty() || allowedTris.find(tCndN) != allowedTris.end()))
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
		std::cout << "face: ";
		for(auto i: f) {
			std::cout << " " << i;
		}
		std::cout <<std::endl;
		visitedTris.insert(triIdx);
		auto fIdx = addFace(f);

		for(auto tidx : trisInFace) {
			auto it = triToFaceMap.find(tidx);
			if(it != triToFaceMap.end()) {
				triToFaceMap.erase(it);
			}
			triToFaceMap.insert({{tidx, fIdx}});
		}
		auto it = faceToTriMap.find(fIdx);
		if(it != faceToTriMap.end()) {
			faceToTriMap.erase(it);
		}
		faceToTriMap.insert({{fIdx, trisInFace}});
	}
}

bool Broker::addTriToFace(long int tidx, Face& f) {
	auto t = tri.getTriangle(tidx);
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
		if(f.size() < 2) {continue;}
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
