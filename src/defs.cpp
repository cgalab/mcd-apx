#include "defs.h"

std::ostream& operator<< (std::ostream& os, const pnt& p) {
	os <<   "x: "  << p.x
	   << ", y: "  << p.y
	   << ", id: " << p.id
	   << ", in: " << std::boolalpha << p.in;
	return os;
}

