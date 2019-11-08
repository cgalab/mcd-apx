#pragma once

#include <functional>
#include <queue>
#include <array>
#include <vector>
#include <string>
#include <assert.h>
#include <iostream>
#include <random>
#include <set>

#include <string.h>

#include "defs.h"
#include "numerics.h"
#include "data.h"

/* interact with triangle */
#ifdef __cplusplus
extern "C" {
	#include "../triangle/triangle.h"
}
#endif

using Edge = std::array<pnt, 2>;

class Triangle {
public:
	Triangle(int id=0, int a=0, int b=0, int c=0, int nAB=0, int nBC=0, int nCA=0)
	: id(id)
	, a(a), b(b), c(c)
	, nAB(nAB), nBC(nBC), nCA(nCA)	{}
	~Triangle() {}

	void setClassic(int id_, int a_, int b_, int c_, int nBC_, int nCA_, int nAB_) {
		id=id_;a=a_;b=b_;c=c_;nAB=nAB_;nBC=nBC_;nCA=nCA_;
	}

	void setNeighbors(int nBC_, int nCA_, int nAB_) {nAB=nAB_;nBC=nBC_;nCA=nCA_;}

	long id;
	long a, b, c;
	long nAB, nBC, nCA;

	int diagonalNeighbor(int i) const {
		if(i == a) { return nBC; }
		if(i == b) { return nCA; }
		if(i == c) { return nAB; }
		return NIL;
	}

	long getCWCorner(long i) const {
		assert(i!=a && i!=b && i!=c);
		if(i == a) {return c;}
		if(i == b) {return a;}
		if(i == c) {return b;}
		return MAX;
	}
	long getCCWCorner(int i) const {
		assert(i!=a && i!=b && i!=c);
		if(i == a) {return b;}
		if(i == b) {return c;}
		if(i == c) {return a;}
		return MAX;
	}

	long getThirdIndex(long int i, long int j) const {
		for(auto idx : {a,b,c}) {
			if(idx != i && idx != j) {
				return idx;
			}
		}
		return NIL;
	}

	bool hasIndex(long int i) const {return i == a || i == b || i == c;}

	void updateNeibhorFromTo(int nOld, int nNew) {
		if(nAB == nOld) {nAB=nNew;}
		else if(nBC == nOld) {nBC=nNew;}
		else if(nCA == nOld) {nCA=nNew;}
	}

	bool isValid() { return a != b && b != c && a != c; }

	void printObjFace() {std::cout << "f " << a+1 << " " << b+1 << " " << c+1 << std::endl;}

	friend std::ostream& operator<<(std::ostream& os, const Triangle& dt);
	friend bool operator==(const Triangle& a, const Triangle& b);
	friend bool operator!=(const Triangle& a, const Triangle& b);
};

using Triangles = std::vector<Triangle>;

class Tri {
public:
	Tri() {
		/* z..index vertices from zero
		 * p..Triangulate (PSLG)
		 * n..compute triangle neighbors
		 * V..Verbose
		 * Q..quiet
		 * */
		std::string mysw = "pznYQce";
		triswitches = new char[mysw.length()+1];
		strcpy(triswitches,mysw.c_str());
	}

	~Tri() {delete triswitches; }

	void runTriangle(pnt* pnts, int num_pnts, Pnts holePnts, const Onions& onions);

	void resetForSortedFlipping();

	bool isEmpty() { return tOUT.numberofpoints == 0; }
	bool isTriangulationDone() { return triangulationDone; }

	Triangle getTriangle(int idx) const {
		assert(idx < (long)tOUT.numberoftriangles);
		return Triangle(idx,tOUT.trianglelist[idx*3],tOUT.trianglelist[idx*3 + 1],tOUT.trianglelist[idx*3 + 2],
					    tOUT.neighborlist[idx*3 + 2],tOUT.neighborlist[idx*3],tOUT.neighborlist[idx*3 + 1]);
	}

	Triangle findTriangleWithCorner(const int idx) const;

	void writeBack(const Triangle& tri);

	Triangle getNextCCWTriangleAroundVertex(const Triangle& tri, long vertex) const;
	Triangle getNextCWTriangleAroundVertex(const Triangle& tri, long vertex) const;


	bool hasCorner(const Triangle& tri, const int a) const {
		return tri.a == a || tri.b == a || tri.c == a;
	}

	bool isOnVertices(const Triangle& tri, const int a, const int b, const int c) const {
		return hasCorner(tri,a) &&  hasCorner(tri,b) && hasCorner(tri,c);
	}


	std::array<Edge,3> getTriangleEdges(int idx) const {
		Triangle t = getTriangle(idx);
		Edge a{{P(t.a), P(t.b)}};
		Edge b{{P(t.b), P(t.c)}};
		Edge c{{P(t.c), P(t.a)}};
		return {{a,b,c}};
	}

	Edge getEdge(int idx) const {
		int PaIdx = tOUT.edgelist[idx*2];
		int PbIdx = tOUT.edgelist[idx*2+1];

		pnt a{tOUT.pointlist[PaIdx], tOUT.pointlist[PaIdx+1], PaIdx, false};
		pnt b{tOUT.pointlist[PbIdx], tOUT.pointlist[PbIdx+1], PbIdx, false};
		return Edge{{a,b}};
	}

	pnt P(int idx) const {
		return pnt{tOUT.pointlist[idx*2], tOUT.pointlist[idx*2 + 1],idx*2,false};
	}

	const triangulateio* getTriangleData() const { return &tOUT; }

	void repairTriangulationOn(std::list<long> tris, const int vertex);

	void reflexSensitiveFlipping(Triangle tri, const int vertex);

	void flipPair(Triangle& ta, Triangle& tb);
	void flipPair(int a, int b) {flipPair(triangles[a],triangles[b]);}

	std::array<long,2> getCommonPair(const Triangle& ta, const Triangle& tb) const;

	int getMissingCorner(const Triangle& tri, int x, int y) const {
		std::array<long,3> array = {{tri.a,tri.b,tri.c}};
		for(auto i : array) {
			if(i != x && i != y) {
				return i;
			}
		}
		assert(false);
		return 0;
	}

	Triangle getNeighborWithoutVertex(const Triangle t, const int vertex);

	bool isConvexQuad(const Triangle& ta, const Triangle& tb) const;

	void setConfig(rt_options* rt_opt_)  {rt_opt = rt_opt_;}

	void printTriangles() const;

	Triangles triangles;

private:
	void filltriangulateioIn(pnt* pnts, int num_pnts, Pnts holePnts, const Onions& onions, triangulateio& tri);
	void inittriangulateioOut(int num_pnts, const Onions& onions, triangulateio& tri);

	triangulateio triangleIN, tOUT, vorout;
	rt_options* rt_opt = nullptr;

	char *triswitches;
	bool triangulationDone = false;

	bool randomSelection 			= false;
	std::random_device rd;
};
