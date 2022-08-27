#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include <cmath>
#include"functions.h"
#include"class.h"

Node* add_nodes(int name, std::vector<Node*>& nodes) {
	if (nodes.size() == 0) {
		Node* node;
		node = new Node(name);
		nodes.push_back(node);
		return node;
	}
	else {
		bool repeat = false;
		int number = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if (name == nodes[i]->get_node_orig()) {
				repeat = true;
				number = i;
			}
		}
		if (repeat == true) {
			return nodes[number];
		}
		else {
			Node* node;
			node = new Node(name);
			nodes.push_back(node);
			return node;
		}
	}
}


