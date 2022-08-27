#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include <cmath>
#include"class.h"
void open(std::vector< Current_source*>& Isourcedc, std::vector<Capacitor*>capacitors, std::vector<int> allnodes) {
	for (int i = 0; i < Isourcedc.size(); i++) {
		int node1 = Isourcedc[i]->get_n_pos()->get_node_dc();
		int node2 = Isourcedc[i]->get_n_neg()->get_node_dc();
		for (int c = 0; c < capacitors.size(); c++) {
			int nodec1 = capacitors[c]->get_node1()->get_node_dc();
			int nodec2 = capacitors[c]->get_node2()->get_node_dc();
			std::cout << nodec1 << " " << nodec2 << " ";
			if ((node1 == nodec1) || (node1 == nodec2)) {
				int count = 0;
				for (int x = 0; x < allnodes.size(); x++) {
					if (allnodes[x] == node1) {
						count++;
					}
				}
				if (count < 2) {
					Isourcedc[i]->get_I() = 0;
				}
			}
			else if ((node2 == nodec1) || (node2 == nodec2)) {
				int count = 0;

				for (int x = 0; x < allnodes.size(); x++) {
					if (allnodes[x] == node2) {
						count++;
					}
				}
				if (count < 2) {
					Isourcedc[i]->change_I(0);
					std::cout << Isourcedc[i]->get_I();
				}

			}
		}
	}
}

