
#include "Tri.h"

std::ostream& operator<<(std::ostream& os, const Triangle& dt) {
    os << "(" << dt.id << ") " << dt.a << ',' << dt.b << ',' << dt.c << " N[" << dt.nBC << "," << dt.nCA << "," << dt.nAB << "]";
    return os;
}
bool operator==(const Triangle& a, const Triangle& b) {return a.id == b.id;}
bool operator!=(const Triangle& a, const Triangle& b) {return a.id != b.id;}


void Tri::runTriangle(Data& data) {
	this->data = &data;

	/* init/fill the triangle structure with data */
	filltriangulateioIn( data,triangleIN);
	inittriangulateioOut(data,tOUT 	    );

	/* triangulate */
	triangulate(triswitches,&triangleIN,&tOUT,&vorout);
	triangulationDone = true;
}



void Tri::flipPair(Triangle ta, Triangle tb) {
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
	return CCW(&P(taM),&P(cp[0]),&P(tbM)) && CCW(&P(tbM),&P(cp[1]),&P(taM));
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

	if(!CCW(&P(vect[0]),&P(vect[1]),&P(getMissingCorner(ta,vect[0],vect[1])) ) ) {
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

//bool Tri::isTriOnBoundaryABAndReflexVertexC(const Triangle& tri) const {
//	if((data->hasEdge(tri.a,tri.b) && data->isReflexVertex(tri.c)) ||
//	   (data->hasEdge(tri.b,tri.c) && data->isReflexVertex(tri.a)) ||
//	   (data->hasEdge(tri.c,tri.a) && data->isReflexVertex(tri.b))
//	) {
//		return true;
//	}
//	return false;
//}
//
//bool Tri::isTriIcidentOnReflexVertexAndBoundary(const Triangle& tri) const {
//	if((data->hasEdge(tri.a,tri.b) && (data->isReflexVertex(tri.a) || data->isReflexVertex(tri.b) ) ) ||
//  	   (data->hasEdge(tri.b,tri.c) && (data->isReflexVertex(tri.b) || data->isReflexVertex(tri.c) ) ) ||
//  	   (data->hasEdge(tri.c,tri.a) && (data->isReflexVertex(tri.c) || data->isReflexVertex(tri.a) ) )
//	) {
//		return true;
//	}
//	return false;
//}

void Tri::filltriangulateioIn(Data& data, triangulateio& tri) {
	auto& vertices 		= data.pnts;
	auto  num_vertices  = data.num_pnts;
	auto& polygon  		= data.onionZero;

	tri.numberofpoints	= num_vertices;
	tri.pointlist    	= new double[tri.numberofpoints * 2];
	tri.pointmarkerlist = new int[tri.numberofpoints];

	for(long i = 0; i < num_vertices; ++i) {
		tri.pointlist[2*i]     = vertices[i].x();
		tri.pointlist[2*i + 1] = vertices[i].y();
		/* one is default for CH boundary vertices */
		tri.pointmarkerlist[i] = i+2;
	}

	tri.numberofholes 			= 0;
	tri.holelist				= NULL;
	tri.numberofregions			= 0;
	tri.regionlist				= NULL;
	tri.numberofcorners			= 3;

	tri.numberofpointattributes = 0;
	tri.pointattributelist 		= NULL;

	tri.numberofsegments = polygon.size();
	tri.segmentlist = new int[tri.numberofsegments * 2];
	tri.segmentmarkerlist = new int[tri.numberofsegments];
	for(long i = 0; i < polygon.size(); ++i) {
		tri.segmentlist[2*i]     = polygon[i];
		tri.segmentlist[2*i + 1] = polygon[(i+1)%polygon.size()];

		/* one is default for CH boundary segments */
		tri.segmentmarkerlist[i] = i+2;
	}
}

void Tri::inittriangulateioOut(Data& data, triangulateio& tri) {
	auto& vertices   	  = data.pnts;
	tri.numberofpoints	  = data.num_pnts;
	tri.pointlist    	  = NULL; //new double[tri.numberofpoints * 2];
	tri.pointmarkerlist   = NULL; //new int[tri.numberofpoints];

	assert(data.num_pnts > 2);
	tri.numberoftriangles = data.num_pnts - 2;
	tri.trianglelist 	  = NULL; //new int[tri.numberoftriangles * 3];
	tri.neighborlist 	  = NULL; //new int[tri.numberoftriangles * 3];
	tri.numberofcorners	  = 3;

	tri.numberofedges 	  = 2 * data.num_pnts - 2;
	tri.edgelist 		  = NULL; //new int[tri.numberofedges * 2];
	tri.edgemarkerlist 	  = NULL; //new int[tri.numberofedges];

	auto& polygon  	 	  = data.onionZero;
	tri.numberofsegments  = polygon.size();
	tri.segmentlist 	  = NULL; //new int[tri.numberofsegments * 2];
	tri.segmentmarkerlist = NULL; //new int[tri.numberofsegments];
}


void Tri::printTriangles() const {
	for (long i = 0; i < (long)tOUT.numberoftriangles; ++i) {
		Triangle t = getTriangle(i);
		std::cout << t << " ";
	}
	std::cout << std::endl;
}


void Tri::printEdges() const {
	std::cout << "edgelist idx: ";
	for (long i = 0; i < (long)tOUT.numberofedges*2; ++i) {
		std::cout << tOUT.edgelist[i] << ", ";
	}
	std::cout << std::endl;

	std::cout << std::endl;
	for (long i = 0; i < (long)tOUT.numberofedges; ++i) {
		Edge e = getEdge(i);
		std::cout << e << " / ";
	}
	std::cout << std::endl;
}


