#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include "../eigen/Eigen/Dense"
#include"functions.h"
#include"class.h"
#include"short_circuit.h"
#include"construct_data_from_file.h"
#include"matrix_non_lin.h"
#include"to_equiv.h"
#include"output.h"
#include"open.h"
#include"nodes.h"

int main() {

    std::ifstream file("test/TEST5.txt");
    std::vector<std::string> s;
    std::string str;
    if (!file) {
        std::cout << "Couldn't open the file" << std::endl;
    }
    else {
        std::vector<ConductanceDevice*> devices;
        std::vector<ConductanceDevice*> resistors;
        std::vector<Inductor*> inductors;
        std::vector<Capacitor*> capacitors;
        std::vector<BJT*> BJTs;
        std::vector<Mosfet*> Mosfets;
        std::vector<Diode*> Diodes;
        std::vector<Voltage_source*> Vsourcesdc;
        std::vector<Voltage_source*> Vsourcesac;
        std::vector<Current_source*> Isourcesdc;
        std::vector<Current_source*> Isourcesac;
        std::vector<Node*>nodes;
        std::vector<int>allnodes;
        double fs = 0, fe = 0;
        int interval = 0;

        while (std::getline(file, str)) {

            std::stringstream ss(str);
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> vstrings(begin, end);
            std::vector<std::string> v = vstrings;
            construct_data(fs, fe, interval, inductors, resistors, v, devices, BJTs, Mosfets, Diodes, Vsourcesdc, Vsourcesac, Isourcesdc, Isourcesac, nodes, capacitors, allnodes);
        }
        
        //vector stores all frequency
        std::vector<double> frequencyv;
        to_f(frequencyv, fs, fe, interval);
        //vector stores all voltages
        std::vector<std::complex<double>> voltages;


        // For DC analysis, short AC voltage sources and inductors

        short_v_source_dc(Vsourcesac, nodes);
        short_inductor(inductors, nodes);
        open(Isourcesdc, capacitors, allnodes);

        // For AC analysis, short DC voltage source
        short_v_source_ac(Vsourcesdc, nodes);
        


        // DC ANALYSIS
        std::vector <ConductanceDevice*> large_sig_devices;
        std::vector<Vcontrolled_Isource*> G;
        std::vector<iSource_nonlin*> isource_non_lin;
        std::vector<Resistor_non_lin_equiv*> non_lin_resis_equiv;

        dc_convert(large_sig_devices, G, isource_non_lin, non_lin_resis_equiv, resistors, BJTs, Mosfets, Diodes);

        //set constant for Newton Raphson
        double Eabs = 0.01;
        double Erel = 0.01;


        solve_initial_solution(BJTs, Mosfets, Diodes, nodes);
        for (int i = 0; i < nodes.size(); i++) {
            std::cout << "N00" << nodes[i]->get_node_orig() << ": " << nodes[i]->get_voltage_dc() << std::endl;
        }
        dc_analysis_newton_raphson(Eabs, Erel, nodes, large_sig_devices, G, isource_non_lin, non_lin_resis_equiv, Vsourcesdc, Isourcesdc);


        //enter output node 
        std::string n_ref, n_get;
        std::cout << "please enter the reference node: " << std::endl;
        std::cin >> n_ref;
        std::cout << "please enter the output node: " << std::endl;
        std::cin >> n_get;
        int ref = stoi(n_ref);
        int get = stoi(n_get);

        std::vector<double> magnitudes;
        std::vector<double> phases;

        // ac analysis
        std::vector<ConductanceDevice*> ssem_devices;
        ac_convert(ssem_devices, non_lin_resis_equiv, devices);


        for (int i = 0; i < frequencyv.size(); i++) {
            double w = 2 * M_PI * frequencyv[i];
            ac_analysis_ssem(w, nodes, ssem_devices, non_lin_resis_equiv, G, Vsourcesac, Isourcesac);
            double mag_ref = 0, mag_get = 0, phase_ref = 0, phase_get = 0;
            for (int r = 0; r < nodes.size(); r++) {
                if (nodes[r]->get_node_orig() == ref) {
                    std::complex<double> voltage = nodes[r]->get_voltage_ac();
                    mag_ref = std::abs(voltage);
                    phase_ref = std::arg(voltage);
                }
                else if (nodes[r]->get_node_orig() == get) {
                    std::complex<double> voltage = nodes[r]->get_voltage_ac();
                    mag_get = std::abs(voltage);
                    phase_get = std::arg(voltage);
                }
            }
            double magnitude = mag_get / mag_ref;
            magnitudes.push_back(magnitude);
            double phase = (phase_get - phase_ref) * (360 / (2 * M_PI));
            std::cout << "phaseout: " << phase_get << std::endl;
            std::cout << "phaseref: " << phase_ref << std::endl;
            phases.push_back(phase);
        }

        //output
        std::ofstream outfile;
        outfile.open("simdata.txt");
        outfile.clear();


        if (!outfile.is_open()) {
            std::cout << "error opening file" << std::endl;
            return EXIT_FAILURE;
        }
        outfile << "frequency" <<"\t"<< "magnitude" << "\t" << "phase" << std::endl;
        for (int i = 0; i < frequencyv.size(); i++) {
            outfile << frequencyv[i] << "\t" << magnitudes[i] << "\t" << phases[i] << std::endl;
        }

        outfile.close();
    }
    return 0;
}






