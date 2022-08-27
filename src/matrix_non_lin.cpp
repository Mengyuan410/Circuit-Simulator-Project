#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <complex>
#include "class.h"
#include "../eigen/Eigen/Dense"
#include <algorithm>

//write dc voltage to node objects
void change_dc_voltage_node(std::vector<Node*>& nodes, const Eigen::VectorXcd Vx, const std::vector<int> large_sig_nodes) {
    for (int r = 0; r < nodes.size(); r++) {
        for (int i = 0; i < large_sig_nodes.size(); i++) {
            if (nodes[r]->get_node_dc() == large_sig_nodes[i]) {
                nodes[r]->change_voltage_dc(Vx[i]);
            }
        }
    }
}

//build vectors of different nodes for dc analysis
std::vector<int> build_dc_nodes_vector(const std::vector<ConductanceDevice*> devices, const std::vector<Voltage_source*> Vsourcesdc) {
    std::vector<int> dc_nodes;
    int n = 0;
    while (devices[n]->get_node1()->get_node_dc() == 0) {
        n = n + 1;
    }
    int x = devices[n]->get_node1()->get_node_dc();
    dc_nodes.push_back(x);
    for (int i = 0; i < devices.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < dc_nodes.size() && newnode == true) {
            if (devices[i]->get_node1()->get_node_dc() == dc_nodes[j] || devices[i]->get_node1()->get_node_dc() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            dc_nodes.push_back(devices[i]->get_node1()->get_node_dc());
        }
    }
    for (int i = 0; i < devices.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < dc_nodes.size() && newnode == true) {
            if (devices[i]->get_node2()->get_node_dc() == dc_nodes[j] || devices[i]->get_node2()->get_node_dc() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            dc_nodes.push_back(devices[i]->get_node2()->get_node_dc());
        }
    }
    for (int i = 0; i < Vsourcesdc.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < dc_nodes.size() && newnode == true) {
            if (Vsourcesdc[i]->get_n_pos()->get_node_dc() == dc_nodes[j] || Vsourcesdc[i]->get_n_pos()->get_node_dc() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            dc_nodes.push_back(Vsourcesdc[i]->get_n_pos()->get_node_dc());
        }
    }
    for (int i = 0; i < Vsourcesdc.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < dc_nodes.size() && newnode == true) {
            if (Vsourcesdc[i]->get_n_neg()->get_node_dc() == dc_nodes[j] || Vsourcesdc[i]->get_n_neg()->get_node_dc() == 0) {
                newnode = false;
            }
            j++;
        }

        if (newnode == true) {
            dc_nodes.push_back(Vsourcesdc[i]->get_n_neg()->get_node_dc());
        }
    }
    return dc_nodes;
}

std::vector<int> build_ac_nodes_vector(const std::vector<ConductanceDevice*> devices, const std::vector<Voltage_source*> Vsourcesac) {
    std::vector<int> ac_nodes;
    int n = 0;
    while (devices[n]->get_node1()->get_node_ac() == 0) {
        n = n + 1;
    }
    int x = devices[n]->get_node1()->get_node_ac();
    ac_nodes.push_back(x);
    for (int i = 0; i < devices.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < ac_nodes.size() && newnode == true) {
            if (devices[i]->get_node1()->get_node_ac() == ac_nodes[j] || devices[i]->get_node1()->get_node_ac() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            ac_nodes.push_back(devices[i]->get_node1()->get_node_ac());
        }
    }
    for (int i = 0; i < devices.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < ac_nodes.size() && newnode == true) {
            if (devices[i]->get_node2()->get_node_ac() == ac_nodes[j] || devices[i]->get_node2()->get_node_ac() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            ac_nodes.push_back(devices[i]->get_node2()->get_node_ac());
        }
    }
    for (int i = 0; i < Vsourcesac.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < ac_nodes.size() && newnode == true) {
            if (Vsourcesac[i]->get_n_pos()->get_node_ac() == ac_nodes[j] || Vsourcesac[i]->get_n_pos()->get_node_ac() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            ac_nodes.push_back(Vsourcesac[i]->get_n_pos()->get_node_ac());
        }
    }
    for (int i = 0; i < Vsourcesac.size(); i++) {
        bool newnode = true;
        int j = 0;
        while (j < ac_nodes.size() && newnode == true) {
            if (Vsourcesac[i]->get_n_neg()->get_node_ac() == ac_nodes[j] || Vsourcesac[i]->get_n_neg()->get_node_ac() == 0) {
                newnode = false;
            }
            j++;
        }
        if (newnode == true) {
            ac_nodes.push_back(Vsourcesac[i]->get_n_neg()->get_node_ac());
        }
    }
    return ac_nodes;
}



///////////////////////////////////INITIAL Solution///////////////////////////////////////////////////
void  solve_initial_solution(const std::vector<BJT*> bjt, const std::vector<Mosfet*> mosfets, const std::vector<Diode*> diodes, std::vector<Node*>& nodes) {
    for (int i = 0; i < bjt.size(); i++) {
        if (bjt[i]->get_model() == "NPN") {
            for (int j = 0; j < nodes.size(); j++) {
                if ((bjt[i]->get_base()->get_node_orig() == nodes[j]->get_node_orig()) && (bjt[i]->get_base()->get_node_dc() != 0)) {
                    nodes[j]->change_voltage_dc(0.7);
                }
                if ((bjt[i]->get_collector()->get_node_orig() == nodes[j]->get_node_orig()) && (bjt[i]->get_collector()->get_node_dc() != 0)) {
                    nodes[j]->change_voltage_dc(0.7);
                }
                if ((bjt[i]->get_emitter()->get_node_orig() == nodes[j]->get_node_orig()) && (bjt[i]->get_emitter()->get_node_dc() != 0)) {
                    nodes[j]->change_voltage_dc(0.7);
                }
            }
        }
        else if (bjt[i]->get_model() == "PNP") {
            for (int j = 0; j < nodes.size(); j++) {
                if ((bjt[i]->get_emitter()->get_node_orig() == nodes[j]->get_node_orig()) && (bjt[i]->get_emitter()->get_node_dc() != 0)) {
                    nodes[j]->change_voltage_dc(0.9);
                }
                if ((bjt[i]->get_base()->get_node_orig() == nodes[j]->get_node_orig()) && (bjt[i]->get_base()->get_node_dc() != 0)) {
                    nodes[j]->change_voltage_dc(0.2);
                }
            }
        }
    }
    for (int i = 0; i < diodes.size(); i++) {
        for (int j = 0; j < nodes.size(); j++) {
            if ((diodes[i]->get_anode()->get_node_orig() == nodes[j]->get_node_orig())&&(diodes[i]->get_anode()->get_node_dc() != 0)){
                nodes[j]->change_voltage_dc(0.7+(diodes[i]->get_cathode()->get_voltage_dc()));
            }
        }
    }
    for (int i = 0; i < mosfets.size(); i++) {
        for (int j = 0; j < nodes.size(); j++) {
            if ((mosfets[i]->get_gate()->get_node_orig() == nodes[j]->get_node_orig())&& (mosfets[i]->get_gate()->get_node_orig()!=0)) {
                nodes[j]->change_voltage_dc(mosfets[i]->get_Vt_mosfet() + (mosfets[i]->get_source()->get_voltage_dc()));
            }
        }
    }
}

//////////////////////////////////DC ANALYSIS WITH NON_LINEAR COMPONENTS////////////////////////////////
//build Jacobian matrix
Eigen::MatrixXcd build_jacobian_matrix(const int size, const int gsize, const std::vector<int> large_sig_nodes,
    std::vector<ConductanceDevice*>& devices_large_sig, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv,
    std::vector<Vcontrolled_Isource*>& G, const std::vector<Voltage_source*> Vsourcesdc) {

    //calculate new resisitance from Vx, new voltage stored in nodes
    for (int r = 0; r < non_lin_resis_equiv.size(); r++) {
        std::complex<double> V_a_b_g = 0, V_0_e_s = 0, V_c_c_d = 0;
        V_a_b_g = non_lin_resis_equiv[r]->get_n_a_b_g()->get_voltage_dc();
        if (non_lin_resis_equiv[r]->get_n_0_e_s() != NULL) {
            V_0_e_s = non_lin_resis_equiv[r]->get_n_0_e_s()->get_voltage_dc();
        }
        else {
            V_0_e_s = 0;
        }
        V_c_c_d = non_lin_resis_equiv[r]->get_n_c_c_d()->get_voltage_dc();
        non_lin_resis_equiv[r]->find_conductance(V_a_b_g, V_0_e_s, V_c_c_d);
    }
    //calculate new gm from Vx, new voltage sored in nodes
    for (int r = 0; r < G.size(); r++) {
        std::complex<double> V_a_b_g = 0, V_0_e_s = 0, V_c_c_d = 0;
        V_a_b_g = G[r]->get_n_a_b_g()->get_voltage_dc();
        V_0_e_s = G[r]->get_n_0_e_s()->get_voltage_dc();
        V_c_c_d = G[r]->get_n_c_c_d()->get_voltage_dc();
        G[r]->find_gm(V_a_b_g, V_0_e_s, V_c_c_d);
    }
    Eigen::MatrixXcd jacobian_matrix = Eigen::MatrixXcd::Zero(size, size);
    //build G
    for (int i = 0; i < gsize; i++) {
        for (int j = 0; j < gsize; j++) {
            for (int r = 0; r < devices_large_sig.size(); r++) {
                if ((devices_large_sig[r]->get_node1()->get_node_dc() == large_sig_nodes[i] && devices_large_sig[r]->get_node2()->get_node_dc() == large_sig_nodes[j]) || (devices_large_sig[r]->get_node1()->get_node_dc() == large_sig_nodes[j] && devices_large_sig[r]->get_node2()->get_node_dc() == large_sig_nodes[i])) {
                    jacobian_matrix(i, j) = jacobian_matrix(i, j) - devices_large_sig[r]->get_conductance(0);
                }
            }
        }
    }
    for (int i = 0; i < gsize; i++) {
        for (int r = 0; r < devices_large_sig.size(); r++) {
            if (devices_large_sig[r]->get_node1()->get_node_dc() == large_sig_nodes[i] || devices_large_sig[r]->get_node2()->get_node_dc() == large_sig_nodes[i]) {
                jacobian_matrix(i, i) = jacobian_matrix(i, i) + devices_large_sig[r]->get_conductance(0);
            }
        }
    }
    //add G into matrix(voltage controlled current sources)
    for (int r = 0; r < G.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (G[r]->get_pos()->get_node_dc() == large_sig_nodes[i]) {
                for (int j = 0; j < gsize; j++) {
                    if (G[r]->get_con_pos()->get_node_dc() == large_sig_nodes[j]) {
                        jacobian_matrix(i, j) = jacobian_matrix(i, j) + G[r]->get_gm();
                    }
                    else if (G[r]->get_con_neg()->get_node_dc() == large_sig_nodes[j]) {
                        jacobian_matrix(i, j) = jacobian_matrix(i, j) - G[r]->get_gm();
                    }

                }
            }
            else if (G[r]->get_neg()->get_node_dc() == large_sig_nodes[i]) {
                for (int j = 0; j < gsize; j++) {
                    if (G[r]->get_con_pos()->get_node_dc() == large_sig_nodes[j]) {
                        jacobian_matrix(i, j) = jacobian_matrix(i, j) - G[r]->get_gm();
                    }
                    else if (G[r]->get_con_neg()->get_node_dc() == large_sig_nodes[j]) {
                        jacobian_matrix(i, j) = jacobian_matrix(i, j) + G[r]->get_gm();
                    }

                }
            }
        }
    }

    //build B
    for (int r = 0; r < Vsourcesdc.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (Vsourcesdc[r]->get_n_pos()->get_node_dc() == large_sig_nodes[i]) {
                jacobian_matrix(i, r + gsize) = 1;
            }
            else if (Vsourcesdc[r]->get_n_neg()->get_node_dc() == large_sig_nodes[i]) {
                jacobian_matrix(i, r + gsize) = -1;
            }
        }
    }
    //build C
    for (int r = 0; r < Vsourcesdc.size(); r++) {
        for (int j = 0; j < gsize; j++) {
            if (Vsourcesdc[r]->get_n_pos()->get_node_dc() == large_sig_nodes[j]) {
                jacobian_matrix(r + gsize, j) = 1;
            }
            else if (Vsourcesdc[r]->get_n_neg()->get_node_dc() == large_sig_nodes[j]) {
                jacobian_matrix(r + gsize, j) = -1;
            }
        }
    }
    return jacobian_matrix;
}

//build vi column with non-linear devices V-controlled I source
Eigen::VectorXcd build_dc_visource_column(const int size, const int gsize, const std::vector<int> large_sig_nodes,
    std::vector<iSource_nonlin*>& Isource_non_lin, const std::vector<Voltage_source*> Vsourcesdc, const std::vector<Current_source*>Isourcesdc) {

    //calculate new resisitance code here from Vx, new voltage in the node
    for (int r = 0; r < Isource_non_lin.size(); r++) {
        std::complex<double> V_a_b_g = 0, V_0_e_s = 0, V_c_c_d = 0;
        V_a_b_g = Isource_non_lin[r]->get_n_a_b_g()->get_voltage_dc();
        if (Isource_non_lin[r]->get_n_0_e_s() != NULL) {
            V_0_e_s = Isource_non_lin[r]->get_n_0_e_s()->get_voltage_dc();
        }
        else {
            V_0_e_s = 0;
        }
        V_c_c_d = Isource_non_lin[r]->get_n_c_c_d()->get_voltage_dc();
        Isource_non_lin[r]->find_I(V_a_b_g, V_0_e_s, V_c_c_d);
    }


    Eigen::VectorXcd vicolumn = Eigen::VectorXcd::Zero(size);

    //write voltage source
    for (int i = gsize; i < size; i++) {
        vicolumn(i) = Vsourcesdc[i - gsize]->get_V();
    }
    //write dc current source
    for (int i = 0; i < gsize; i++) {
        for (int r = 0; r < Isourcesdc.size(); r++) {
            if (large_sig_nodes[i] == Isourcesdc[r]->get_n_neg()->get_node_dc()) {
                vicolumn(i) = Isourcesdc[r]->get_I() + vicolumn(i);
            }
            else if (large_sig_nodes[i] == Isourcesdc[r]->get_n_pos()->get_node_dc()) {
                vicolumn(i) = -Isourcesdc[r]->get_I() + vicolumn(i);
            }
        }
    }
    //write current sources provided by non-linear devices
    for (int r = 0; r < Isource_non_lin.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (Isource_non_lin[r]->get_pos()->get_node_dc() == large_sig_nodes[i]) {
                vicolumn(i) = -Isource_non_lin[r]->get_I() + vicolumn(i);
            }
            else if (Isource_non_lin[r]->get_neg()->get_node_dc() == large_sig_nodes[i]) {
                vicolumn(i) = Isource_non_lin[r]->get_I() + vicolumn(i);
            }
        }
    }
    return vicolumn;
}

//Function that can decide whether iteration ends
bool whether_iterate(Eigen::VectorXcd Vx, Eigen::VectorXcd Vx_next, double Eabs, double Erel) {
    bool iterate = false;
    for (int i = 0; i < Vx.size(); i++) {
        std::complex<double> difference = std::abs(Vx_next[i] - Vx[i]);
        std::complex<double> barrier = std::abs(Eabs + Erel * Vx[i]);

        if (std::abs(difference) > std::abs(barrier)) {
            iterate = true;
        }
    }

    return iterate;
}



//Newton Raphson method
void dc_analysis_newton_raphson(const double Eabs, const double Erel, std::vector<Node*>& nodes, std::vector<ConductanceDevice*>& devices_large_sig, std::vector<Vcontrolled_Isource*>& G,
    std::vector<iSource_nonlin*>& Isource_non_lin, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv, const std::vector<Voltage_source*> Vsourcesdc,
    const std::vector<Current_source*> Isourcesdc) {
    std::vector<int> large_sig_nodes = build_dc_nodes_vector(devices_large_sig, Vsourcesdc);
    int size = large_sig_nodes.size() + Vsourcesdc.size();
    int gsize = large_sig_nodes.size();

    Eigen::MatrixXcd jacobian_matrix = build_jacobian_matrix(size, gsize, large_sig_nodes, devices_large_sig, non_lin_resis_equiv, G, Vsourcesdc);

    std::cout << "Initial Jacobian Matrix" << std::endl;
    std::cout << jacobian_matrix << std::endl;

    Eigen::VectorXcd vicolumn = build_dc_visource_column(size, gsize, large_sig_nodes, Isource_non_lin, Vsourcesdc, Isourcesdc);
    std::cout << vicolumn << std::endl;

    std::cout << "Initial RHS Column" << std::endl << vicolumn << std::endl;

    Eigen::VectorXcd Vx = jacobian_matrix.colPivHouseholderQr().solve(vicolumn);

    std::cout << "First Vx Solution: Vx(1)" << std::endl;
    std::cout << Vx << std::endl;

    change_dc_voltage_node(nodes, Vx, large_sig_nodes);
    for (int i = 0; i < nodes.size(); i++) {
        std::cout << "N00" << nodes[i]->get_node_orig() << ": " << nodes[i]->get_voltage_dc() << std::endl;
    }
    bool iterate = true;
    while (iterate == true) {
        jacobian_matrix = build_jacobian_matrix(size, gsize, large_sig_nodes, devices_large_sig, non_lin_resis_equiv, G, Vsourcesdc);
        std::cout << "Jacobian Matrix" << std::endl;
        std::cout << jacobian_matrix << std::endl;
        vicolumn = build_dc_visource_column(size, gsize, large_sig_nodes, Isource_non_lin, Vsourcesdc, Isourcesdc);
        std::cout << vicolumn << std::endl;
        std::cout << "RHS Column" << std::endl << vicolumn << std::endl;
        Eigen::VectorXcd Vx_next = jacobian_matrix.colPivHouseholderQr().solve(vicolumn);
        std::cout << "Vx Solution" << std::endl;
        std::cout << Vx_next << std::endl;
        change_dc_voltage_node(nodes, Vx_next, large_sig_nodes);
        for (int i = 0; i < nodes.size(); i++) {
            std::cout << "N00" << nodes[i]->get_node_orig() << ": " << nodes[i]->get_voltage_dc() << std::endl;
        }
        iterate = whether_iterate(Vx, Vx_next, Eabs, Erel);
        Vx = Vx_next;
    }
}


/////////////////////////////////AC ANALYSIS SSEM/////////////////////////////////////
//build ssem conductance matrix
//non-lin-resis-equiv here has changed nodes for ac analysis i.e. short dc
Eigen::MatrixXcd build_ssem_conductance_matrix(const double omega, int size, int gsize, const std::vector<int> ssemnodes, std::vector<ConductanceDevice*>& devices_small_sig,
    std::vector<Vcontrolled_Isource*>& G, const std::vector<Voltage_source*> Vsourcesac, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv) {
    //calculate R from Vx_dc
    for (int r = 0; r < non_lin_resis_equiv.size(); r++) {
        std::complex<double> V_a_b_g = 0, V_0_e_s = 0, V_c_c_d = 0;
        V_a_b_g = non_lin_resis_equiv[r]->get_n_a_b_g()->get_voltage_dc();
        if (non_lin_resis_equiv[r]->get_n_0_e_s() != NULL) {
            V_0_e_s = non_lin_resis_equiv[r]->get_n_0_e_s()->get_voltage_dc();
        }
        else {
            V_0_e_s = 0;
        }
        V_c_c_d = non_lin_resis_equiv[r]->get_n_c_c_d()->get_voltage_dc();
        non_lin_resis_equiv[r]->find_conductance(V_a_b_g, V_0_e_s, V_c_c_d);
    }
    //calculate new gm from Vx
    for (int r = 0; r < G.size(); r++) {
        std::complex<double> V_a_b_g = 0, V_0_e_s = 0, V_c_c_d = 0;
        V_a_b_g = G[r]->get_n_a_b_g()->get_voltage_dc();
        if (G[r]->get_n_0_e_s() != NULL) {
            V_0_e_s = G[r]->get_n_0_e_s()->get_voltage_dc();
        }
        else {
            V_0_e_s = 0;
        }
        V_c_c_d = G[r]->get_n_c_c_d()->get_voltage_dc();
        G[r]->find_gm(V_a_b_g, V_0_e_s, V_c_c_d);
    }
    Eigen::MatrixXcd gmatrix = Eigen::MatrixXcd::Zero(size, size);
    //build G
    for (int i = 0; i < gsize; i++) {
        for (int j = 0; j < gsize; j++) {
            for (int r = 0; r < devices_small_sig.size(); r++) {
                if ((devices_small_sig[r]->get_node1()->get_node_ac() == ssemnodes[i] && devices_small_sig[r]->get_node2()->get_node_ac() == ssemnodes[j]) || (devices_small_sig[r]->get_node1()->get_node_ac() == ssemnodes[j] && devices_small_sig[r]->get_node2()->get_node_ac() == ssemnodes[i])) {
                    gmatrix(i, j) = gmatrix(i, j) - devices_small_sig[r]->get_conductance(omega);
                }
            }
        }
    }
    for (int i = 0; i < gsize; i++) {
        for (int r = 0; r < devices_small_sig.size(); r++) {
            if (devices_small_sig[r]->get_node1()->get_node_ac() == ssemnodes[i] || devices_small_sig[r]->get_node2()->get_node_ac() == ssemnodes[i]) {
                gmatrix(i, i) = gmatrix(i, i) + devices_small_sig[r]->get_conductance(omega);
            }
        }
    }
    //add G into matrix(voltage controlled current sources)
    for (int r = 0; r < G.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (G[r]->get_pos()->get_node_ac() == ssemnodes[i]) {
                for (int j = 0; j < gsize; j++) {
                    if (G[r]->get_con_pos()->get_node_ac() == ssemnodes[j]) {
                        gmatrix(i, j) = gmatrix(i, j) + G[r]->get_gm();
                    }
                    else if (G[r]->get_con_neg()->get_node_ac() == ssemnodes[j]) {
                        gmatrix(i, j) = gmatrix(i, j) - G[r]->get_gm();
                    }

                }
            }
            else if (G[r]->get_neg()->get_node_ac() == ssemnodes[i]) {
                for (int j = 0; j < gsize; j++) {
                    if (G[r]->get_con_pos()->get_node_ac() == ssemnodes[j]) {
                        gmatrix(i, j) = gmatrix(i, j) - G[r]->get_gm();
                    }
                    else if (G[r]->get_con_neg()->get_node_ac() == ssemnodes[j]) {
                        gmatrix(i, j) = gmatrix(i, j) + G[r]->get_gm();
                    }

                }
            }
        }
    }
    //build B
    for (int r = 0; r < Vsourcesac.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (Vsourcesac[r]->get_n_pos()->get_node_ac() == ssemnodes[i]) {
                gmatrix(i, r + gsize) = 1;
            }
            else if (Vsourcesac[r]->get_n_neg()->get_node_ac() == ssemnodes[i]) {
                gmatrix(i, r + gsize) = -1;
            }
        }
    }

    //build C
    for (int r = 0; r < Vsourcesac.size(); r++) {
        for (int j = 0; j < gsize; j++) {
            if (Vsourcesac[r]->get_n_pos()->get_node_ac() == ssemnodes[j]) {
                gmatrix(r + gsize, j) = 1;
            }
            else if (Vsourcesac[r]->get_n_neg()->get_node_ac() == ssemnodes[j]) {
                gmatrix(r + gsize, j) = -1;
            }
        }
    }
    return gmatrix;
}

Eigen::VectorXcd build_ssem_visource_column(const std::vector<int> ssemnodes, const int size, const int gsize,
    const std::vector<Voltage_source*> Vsourcesac, const std::vector<Current_source*>Isourcesac) {
    Eigen::VectorXcd vicolumn = Eigen::VectorXcd::Zero(size);
    //write voltage source
    for (int i = 0; i < Vsourcesac.size(); i++) {
        vicolumn(i + gsize) = Vsourcesac[i]->get_V();
    }
    //write current voltage
    for (int i = 0; i < gsize; i++) {
        for (int r = 0; r < Isourcesac.size(); r++) {
            if (ssemnodes[i] == Isourcesac[r]->get_n_neg()->get_node_ac()) {
                vicolumn(i) = Isourcesac[r]->get_I() + vicolumn(i);
            }
            else if (ssemnodes[i] == Isourcesac[r]->get_n_pos()->get_node_ac()) {
                vicolumn(i) = -Isourcesac[r]->get_I() + vicolumn(i);
            }
        }
    }
    return vicolumn;
}

//write ac voltage to node objects
void change_ac_voltage_node(std::vector<Node*>& nodes, const Eigen::VectorXcd Vx, const std::vector<int> ssemnodes, const int gsize) {
    for (int r = 0; r < nodes.size(); r++) {
        for (int i = 0; i < gsize; i++) {
            if (nodes[r]->get_node_ac() == ssemnodes[i]) {
                nodes[r]->change_voltage_ac(Vx[i]);
            }
        }
    }
}
//solve Ax=b
void ac_analysis_ssem(double omega, std::vector<Node*>& nodes, std::vector<ConductanceDevice*>& devices_small_sig, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv,
    std::vector<Vcontrolled_Isource*>& G, const std::vector<Voltage_source*> Vsourcesac, const std::vector<Current_source*>Isourcesac) {
    std::vector<int> ssemnodes = build_ac_nodes_vector(devices_small_sig, Vsourcesac);
    int size = ssemnodes.size() + Vsourcesac.size();
    int gsize = ssemnodes.size();
    Eigen::MatrixXcd gmatrix = build_ssem_conductance_matrix(omega, size, gsize, ssemnodes, devices_small_sig, G, Vsourcesac, non_lin_resis_equiv);
    std::cout << "gmatrix" << std::endl << gmatrix << std::endl;
    Eigen::VectorXcd vicolumn = build_ssem_visource_column(ssemnodes, size, gsize, Vsourcesac, Isourcesac);
    std::cout << "vicolumn" << std::endl << vicolumn << std::endl;
    Eigen::VectorXcd Vx = gmatrix.colPivHouseholderQr().solve(vicolumn);
    std::cout << "Vx" << std::endl << Vx << std::endl;
    change_ac_voltage_node(nodes, Vx, ssemnodes, gsize);
    for (int i = 0; i < nodes.size(); i++) {
        std::cout << "N00" << nodes[i]->get_node_orig() << ": " << nodes[i]->get_voltage_ac() << std::endl;
    }
}







