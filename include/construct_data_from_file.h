#ifndef CONSTRUCT_H
#define CONSTRUCT_H

#include "class.h"
void construct_data(double& fs, double& fe, int& interval, std::vector<Inductor*>& inductors, std::vector<ConductanceDevice*>& resistors, const std::vector<std::string> v, std::vector<ConductanceDevice*>& devices, std::vector<BJT*>& BJTs,
    std::vector<Mosfet*>& Mosfets, std::vector<Diode*>& Diodes, std::vector<Voltage_source*>& DCsources, std::vector<Voltage_source*>& ACsources,
    std::vector<Current_source*>& Isources, std::vector<Current_source*>& Isourcesac, std::vector<Node*>& nodes, std::vector<Capacitor*>& capacitors, std::vector<int>& allnodes);
#endif

