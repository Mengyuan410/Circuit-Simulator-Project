#include "class.h"

// DC:
// - short AC V source
// - short inductor


// AC:
// - short DC V source



// For AC: loop through 'nodes' to find node connected to the specified old node, and change it to the new node
void convert_node_ac(std::vector<Node*>& nodes, int old_node, int new_node) {
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_node_ac() == old_node) {
            nodes[i]->change_node_ac(new_node);
        }
    }
}

// For DC: loop through 'nodes' to find node connected to the specified old node, and change it to the new node
void convert_node_dc(std::vector<Node*>& nodes, int old_node, int new_node) {
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_node_dc() == old_node) {
            nodes[i]->change_node_dc(new_node);
        }
    }
}


// Convert voltage sources to short circuit, so the 2 nodes of the voltage source will become the same
// For DC: short AC voltage sources
void short_v_source_dc(std::vector<Voltage_source*> vsources, std::vector<Node*>& nodes) {
    int pos_node, neg_node;
    for (int i = 0; i < vsources.size(); i++) {
        pos_node = vsources[i]->get_n_pos()->get_node_orig();
        neg_node = vsources[i]->get_n_neg()->get_node_orig();
        // if one of the voltage source terminals is connect to ground, make all component connect to the source to ground
        if (pos_node == 0) {
            convert_node_dc(nodes, neg_node, 0);
        }
        else if (neg_node == 0) {
            convert_node_dc(nodes, pos_node, 0);
        }
        // if the voltage source is connect to two non-reference nodes, change node2 to node 1
        else {
            convert_node_dc(nodes, neg_node, pos_node);
        }
    }
}

// For DC: short inductors
void short_inductor(std::vector<Inductor*> inductors, std::vector<Node*>& nodes) {
    int node1, node2;
    for (int i = 0; i < inductors.size(); i++) {
        node1 = inductors[i]->get_node1()->get_node_orig();
        node2 = inductors[i]->get_node2()->get_node_orig();
        // if one of the inductor's terminal is connect to ground, make all component connect to the source to ground
        if (node1 == 0) {
            convert_node_dc(nodes, node2, 0);
        }
        else if (node2 == 0) {
            convert_node_dc(nodes, node1, 0);
        }
        // if the inductor is connect to two non-reference nodes, change node2 to node 1
        else {
            convert_node_dc(nodes, node2, node1);
        }
    }
}


// Convert voltage sources to short circuit, so the 2 nodes of the voltage source will become the same
// For AC: short DC voltage sources
void short_v_source_ac(std::vector<Voltage_source*> vsources, std::vector<Node*>& nodes) {
    int pos_node, neg_node;
    for (int i = 0; i < vsources.size(); i++) {
        pos_node = vsources[i]->get_n_pos()->get_node_orig();
        neg_node = vsources[i]->get_n_neg()->get_node_orig();
        // if one of the voltage source terminals is connect to ground, make all component connect to the source to ground
        if (pos_node == 0) {
            convert_node_ac(nodes, neg_node, 0);
        }
        else if (neg_node == 0) {
          
            convert_node_ac(nodes, pos_node, 0);
        }
        // if the voltage source is connect to two non-reference nodes, change node2 to node 1
        else {
            convert_node_ac(nodes, neg_node, pos_node);
        }
    }
}




