//
//  to_small_sig_equiv.h
//  1st_year_project
//
//  Created by  Cathy Ding on 11/05/2021.


#ifndef to_small_sig_equiv_h
#define to_small_sig_equiv_h
#include "class.h"

ConductanceDevice* D_to_R(Diode d, double VD);
std::vector<ConductanceDevice*> Q_to_R(BJT q, double VBE);
ConductanceDevice* M_to_R(Mosfet m, double GS,double VDS);
Vcontrolled_Isource* Q_to_G(BJT q, double Ic);
Vcontrolled_Isource* M_to_G(Mosfet m, double Id);
void convert_to_R_G(std::vector<ConductanceDevice*>& devices_small_sig, std::vector<Vcontrolled_Isource*>& G_small_sig, std::vector<ConductanceDevice*>& devices, std::vector<BJT*>& BJTs, std::vector<Mosfet*>& Mosfets, std::vector<Diode*>& Diodes, std::vector<double>& current);
void change_node(std::vector<ConductanceDevice*>& devices_small_sig, std::vector<Vcontrolled_Isource*>& G_small_sig, std::vector<Voltage_source_AC*>& ACsources, std::vector<Current_source_AC*>& Isourcesac, int old_node, int new_node);
void short_v_source(std::vector<Voltage_source_DC*> DCsources, std::vector<ConductanceDevice*>& devices_small_sig, std::vector<Vcontrolled_Isource*>& G_small_sig, std::vector<Voltage_source_AC*>& ACsources, std::vector<Current_source_AC*>& Isourcesac);
#endif
