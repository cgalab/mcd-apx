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
	for(auto s : sets) {
		s.printLayer(0);
	}
}

void Broker::printSets() const {
	for(auto s : sets) {
		s.printPnts();
	}
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

void Data::printLayer(const int idx) const {
	std::cout << "idx: " << idx << " num_layers: " << num_layers;
	assert(idx < num_layers);
	auto l = layers[idx];
	auto startNode = nodes[l.nde];
	node *nIt = &startNode;
	std::cout << "print layer " << idx << ", size: " << l.num << std::endl;
	int cnt = 0;
	do {
		auto p = pnts[nIt->vtx];
		std::cout << "idx: " << cnt++ << ", id: " << p.id  << " x: " << p.x << " y: " << p.y << std::endl;
		nIt = &nodes[nIt->next];
	} while(nIt != &startNode);
}



