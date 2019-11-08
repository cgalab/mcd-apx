#include <algorithm>    // std::for_each

#include "Tri.h"

std::ostream& operator<<(std::ostream& os, const Triangle& dt) {
    os << "(" << dt.id << ") " << dt.a << ',' << dt.b << ',' << dt.c << " N[" << dt.nBC << "," << dt.nCA << "," << dt.nAB << "]";
    return os;
}
bool operator==(const Triangle& a, const Triangle& b) {return a.id == b.id;}
bool operator!=(const Triangle& a, const Triangle& b) {return a.id != b.id;}


void Tri::runTriangle(pnt* pnts, int num_pnts, Pnts holePnts, const Onions& onions) {
	/* init/fill the triangle structure with data */
	filltriangulateioIn( pnts,num_pnts,holePnts,onions,triangleIN);
	inittriangulateioOut(     num_pnts,onions,tOUT 	    );

	/* triangulate */
	triangulate(triswitches,&triangleIN,&tOUT,&vorout);
	triangulationDone = true;

	for(int i = 0; i < tOUT.numberoftriangles; ++i) {
		triangles.push_back(getTriangle(i));
	}
}



void Tri::flipPair(Triangle& ta, Triangle& tb) {
	auto cp = getCommonPair(ta,tb);
	long taM = getMissingCorner(ta,cp[0],cp[1]);
	long tbM = getMissingCorner(tb,cp[0],cp[1]);

	Triangle taNew(ta.id, taM,cp[0],tbM, NIL,NIL,NIL); //ta.diagonalNeighbor(cp[0]), tb.diagonalNeighbor(cp[0]) ,tb.id );
	Triangle tbNew(tb.id, tbM,cp[1],taM, NIL,NIL,NIL); //tb.diagonalNeighbor(cp[1]), ta.diagonalNeighbor(cp[1]), ta.id );

	taNew.setNeighbors(tb.diagonalNeighbor(cp[1]),tb.id,ta.diagonalNeighbor(cp[1]));
	tbNew.setNeighbors(ta.diagonalNeighbor(cp[0]),ta.id,tb.diagonalNeighbor(cp[0]));

	writeBack(taNew);
	writeBack(tbNew);

	/* repair back links of neighbours */
	std::vector<Triangle> nUpdate;
	for(auto nid : { tb.diagonalNeighbor(cp[1]), ta.diagonalNeighbor(cp[1]) } ) {
		if(nid != NIL) {
			auto tN = getTriangle(nid);
			tN.updateNeibhorFromTo(tb.id,ta.id);
			nUpdate.push_back(tN);
		}
	}
	for(auto nid : { ta.diagonalNeighbor(cp[0]), tb.diagonalNeighbor(cp[0]) } ) {
		if(nid != NIL) {
			auto tN = getTriangle(nid);
			tN.updateNeibhorFromTo(ta.id,tb.id);
			nUpdate.push_back(tN);
		}
	}
	for(auto t : nUpdate) {
		writeBack(t);
	}
}

bool Tri::isConvexQuad(const Triangle& ta, const Triangle& tb) const {
	auto cp = getCommonPair(ta,tb);
	long taM = getMissingCorner(ta,cp[0],cp[1]);
	long tbM = getMissingCorner(tb,cp[0],cp[1]);
	auto Pa = P(taM); auto Pb = P(tbM);
	auto Pc1 = P(cp[0]); auto Pc2 = P(cp[1]);
	return CCW(&Pa,&Pc1,&Pb) && CCW(&Pb,&Pc2,&Pa);
}

void Tri::writeBack(const Triangle& tri) {
	tOUT.trianglelist[tri.id*3    ] = tri.a;
	tOUT.trianglelist[tri.id*3 + 1] = tri.b;
	tOUT.trianglelist[tri.id*3 + 2] = tri.c;
	tOUT.neighborlist[tri.id*3    ] = tri.diagonalNeighbor(tri.a);
	tOUT.neighborlist[tri.id*3 + 1] = tri.diagonalNeighbor(tri.b);
	tOUT.neighborlist[tri.id*3 + 2] = tri.diagonalNeighbor(tri.c);
}

/* returns the common two corners of ta, tb in the CCW order of 'ta' */
std::array<long,2> Tri::getCommonPair(const Triangle& ta, const Triangle& tb) const {
	std::vector<long> vect;
	if(hasCorner(ta,tb.c)) {vect.push_back(tb.c);}
	if(hasCorner(ta,tb.b)) {vect.push_back(tb.b);}
	if(hasCorner(ta,tb.a)) {vect.push_back(tb.a);}
	assert(vect.size() == 2);

	auto Pv0 = P(vect[0]); auto Pv1 = P(vect[1]);
	auto Pv2 = P(getMissingCorner(ta,vect[0],vect[1]));
	if(!CCW(&Pv0,&Pv1,&Pv2) ) {
		std::swap(vect[0],vect[1]);
	}

	return {{ vect[0],vect[1] }};
}

/* CCW !! rotation */
Triangle Tri::getNextCCWTriangleAroundVertex(const Triangle& tri, long vertex) const {
	if(vertex == tri.a && tri.nCA != NIL) {
		return getTriangle(tri.nCA);
	} else if(vertex == tri.b && tri.nAB != NIL) {
		return getTriangle(tri.nAB);
	} else if(vertex == tri.c && tri.nBC != NIL) {
		return getTriangle(tri.nBC);
	}
	assert(false);
	return Triangle();
}

/* CW !! rotation */
Triangle Tri::getNextCWTriangleAroundVertex(const Triangle& tri, long vertex) const {
	if(vertex == tri.a && tri.nAB != NIL) {
		return getTriangle(tri.nAB);
	} else if(vertex == tri.b && tri.nBC != NIL) {
		return getTriangle(tri.nBC);
	} else if(vertex == tri.c && tri.nCA != NIL) {
		return getTriangle(tri.nCA);
	}
	assert(false);
	return Triangle();
}

void Tri::filltriangulateioIn(pnt* pnts, int num_pnts, Pnts holePnts, const Onions& onions, triangulateio& tri) {
	tri.numberofpoints	= num_pnts;
	tri.pointlist    	= new double[tri.numberofpoints * 2];
	tri.pointmarkerlist = new int[tri.numberofpoints];

	int polygonSize = 0;
	std::for_each(onions.begin(),onions.end(),[&](OnionZero o){polygonSize+=o.size();});

	for(long i = 0; i < num_pnts; ++i) {
		tri.pointlist[2*i]     = pnts[i].x;
		tri.pointlist[2*i + 1] = pnts[i].y;
		/* one is default for CH boundary vertices */
		tri.pointmarkerlist[i] = i+2;
	}

	tri.numberofholes 			= holePnts.size();
	tri.holelist				= new double[tri.numberofholes * 2];
	for(long i = 0; i < tri.numberofholes; ++i) {
		tri.holelist[2*i]     = holePnts[i].x;
		tri.holelist[2*i + 1] = holePnts[i].y;
	}

	tri.numberofregions			= 0;
	tri.regionlist				= NULL;
	tri.numberofcorners			= 3;

	tri.numberofpointattributes = 0;
	tri.pointattributelist 		= NULL;

	tri.numberofsegments = polygonSize;
	tri.segmentlist = new int[tri.numberofsegments * 2];
	tri.segmentmarkerlist = new int[tri.numberofsegments];

	int segCnt = 0;
	for(auto& onion : onions) {
		auto it = onion.begin();
		do {
			tri.segmentlist[2*segCnt] = *it;
			++it;
			if(it == onion.end()) {it = onion.begin();}
			tri.segmentlist[2*segCnt + 1] = *it;
			tri.segmentmarkerlist[segCnt] = segCnt+2;

			++segCnt;
		} while(it != onion.begin());
	}
}

void Tri::inittriangulateioOut(int num_pnts, const Onions& onions, triangulateio& tri) {
	tri.numberofpoints	  = num_pnts;
	tri.pointlist    	  = NULL;
	tri.pointmarkerlist   = NULL;

	assert(num_pnts > 2);
	tri.numberoftriangles = num_pnts - 2;
	tri.trianglelist 	  = NULL;
	tri.neighborlist 	  = NULL;
	tri.numberofcorners	  = 3;

	tri.numberofedges 	  = 2 * num_pnts - 2;
	tri.edgelist 		  = NULL;
	tri.edgemarkerlist 	  = NULL;

	int polygonSize = 0;
	std::for_each(onions.begin(),onions.end(),[&](OnionZero o){polygonSize+=o.size();});

	tri.numberofsegments  = polygonSize;
	tri.segmentlist 	  = NULL;
	tri.segmentmarkerlist = NULL;
}


void Tri::printTriangles() const {
	for (long i = 0; i < (long)tOUT.numberoftriangles; ++i) {
		Triangle t = getTriangle(i);
		std::cout << t << " ";
	}
	std::cout << std::endl;
}



