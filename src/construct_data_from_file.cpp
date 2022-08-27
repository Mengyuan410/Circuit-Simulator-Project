#define _USE_MATH_DEFINES

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
#include <complex>
#include"nodes.h"



void construct_data(double& fs, double& fe, int& interval, std::vector<Inductor*>& inductors, std::vector<ConductanceDevice*>& resistors, const std::vector<std::string> v, std::vector<ConductanceDevice*>& devices, std::vector<BJT*>& BJTs,
    std::vector<Mosfet*>& Mosfets, std::vector<Diode*>& Diodes, std::vector<Voltage_source*>& DCsources, std::vector<Voltage_source*>& ACsources,
    std::vector<Current_source*>& Isources, std::vector<Current_source*>& Isourcesac, std::vector<Node*>& nodes, std::vector<Capacitor*>& capacitors, std::vector<int>& allnodes) {
    Resistor* resistor;
    Inductor* inductor;
    Capacitor* capacitor;
    BJT* bjt;
    Diode* diode;
    Mosfet* mosfet;
    Current_source* isource;
    Current_source* isourceac;
    Voltage_source* acsource;
    Voltage_source* dcsource;
    Node* n1;
    Node* n2;
    Node* n3;
    std::string name = v[0];
    std::string type;
    type = name[0];

    if (type == "R") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        double R = get_value(v[3]);
        resistor = new Resistor(R, n1, n2);
        resistors.push_back(resistor);
        devices.push_back(resistor);
        allnodes.push_back(get_nodes(v[1]));
        allnodes.push_back(get_nodes(v[2]));
    }
    else if (type == "C") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        double C = get_value(v[3]);
        capacitor = new Capacitor(C, n1, n2);
        devices.push_back(capacitor);
        capacitors.push_back(capacitor);
    }
    else if (type == "L") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        double L = get_value(v[3]);
        inductor = new Inductor(L, n1, n2);
        inductors.push_back(inductor);
        devices.push_back(inductor);
        allnodes.push_back(get_nodes(v[1]));
        allnodes.push_back(get_nodes(v[2]));
    }
    else if (type == "D") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        diode = new Diode(n1, n2, "N");
        Diodes.push_back(diode);
        allnodes.push_back(get_nodes(v[1]));
        allnodes.push_back(get_nodes(v[2]));

    }
    else if (type == "Q") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        n3 = add_nodes(get_nodes(v[3]), nodes);
        std::string type = v[4];
        bjt = new BJT(n1, n2, n3, type);
        BJTs.push_back(bjt);
        allnodes.push_back(get_nodes(v[1]));
        allnodes.push_back(get_nodes(v[2]));
        allnodes.push_back(get_nodes(v[3]));
    }
    else if (type == "M") {
        n1 = add_nodes(get_nodes(v[1]), nodes);
        n2 = add_nodes(get_nodes(v[2]), nodes);
        n3 = add_nodes(get_nodes(v[3]), nodes);
        std::string type = v[4];
        mosfet = new Mosfet(n1, n2, n3, type);
        Mosfets.push_back(mosfet);
        allnodes.push_back(get_nodes(v[1]));
        allnodes.push_back(get_nodes(v[2]));
        allnodes.push_back(get_nodes(v[3]));
    }
    else if (name == ".ac") {
        fs = get_value(v[3]);
        fe = get_value(v[4]);
        interval = get_value(v[2]);
    }
    else if (type == "V") {
        if (v[3] == "AC") {
            n1 = add_nodes(get_nodes(v[1]), nodes);
            n2 = add_nodes(get_nodes(v[2]), nodes);
            double amplitude = get_value(v[4]);
            double phase = stod(v[5]);
            std::complex<double> V;
            V.real(amplitude * cos(phase * M_PI / 180));
            V.imag(amplitude * sin(phase * M_PI / 180));
            acsource = new Voltage_source(n1, n2, V);
            ACsources.push_back(acsource);
            allnodes.push_back(get_nodes(v[1]));
            allnodes.push_back(get_nodes(v[2]));

        }
        else {
            n1 = add_nodes(get_nodes(v[1]), nodes);
            n2 = add_nodes(get_nodes(v[2]), nodes);
            double volt = get_value(v[3]);
            dcsource = new Voltage_source(n1, n2, volt);
            DCsources.push_back(dcsource);
            allnodes.push_back(get_nodes(v[1]));
            allnodes.push_back(get_nodes(v[2]));
        }
    }
    else if (type == "I") {
        if (v[3] == "AC") {
            n1 = add_nodes(get_nodes(v[1]), nodes);
            n2 = add_nodes(get_nodes(v[2]), nodes);
            double amplitude = get_value(v[4]);
            double phase = stod(v[5]);
            std::complex<double> I;
            I.real(amplitude * cos(phase * M_PI / 180));
            I.imag(amplitude * sin(phase * M_PI / 180));
            isourceac = new Current_source(n1, n2, I);
            Isourcesac.push_back(isourceac);
            allnodes.push_back(get_nodes(v[1]));
            allnodes.push_back(get_nodes(v[2]));
        }
        else {
            n1 = add_nodes(get_nodes(v[1]), nodes);
            n2 = add_nodes(get_nodes(v[2]), nodes);
            double current = get_value(v[3]);
            isource = new Current_source(n1, n2, current);
            Isources.push_back(isource);
            allnodes.push_back(get_nodes(v[1]));
            allnodes.push_back(get_nodes(v[2]));
        }
    }


}



