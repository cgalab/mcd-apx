#include "data.h"

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
	std::cout << "merge: " << std::endl;
	for(auto s : sets) {
		s.printLayer();
	}
}

void Broker::printSets() const {
	for(auto s : sets) {
		s.printPnts();
	}
}


void Data::backupOnionZoro(int idx) {
	auto l = layers[idx];
	node startNode = nodes[l.nde];
	node nIt = startNode;
	int cnt = 0;
	do {
		std::cout << nIt.vtx << " ";
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

void Data::printLayer() const {
	std::cout << "size of onion: " << onionZero.size() << std::endl;
	for(auto n : onionZero) {
		std::cout << "id: " << n.vtx  << std::endl;
	}
}



