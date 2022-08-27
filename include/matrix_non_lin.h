
#ifndef matrix_non_lin_h
#define matrix_non_lin_h
#include "class.h"
#include "../eigen/Eigen/Dense"

void  solve_initial_solution(const std::vector<BJT*> bjt, const std::vector<Mosfet*> mosfets, const std::vector<Diode*> diodes, std::vector<Node*>& nodes);
void dc_analysis_newton_raphson(const double Eabs, const double Erel, std::vector<Node*>& nodes, std::vector<ConductanceDevice*>& devices_large_sig, std::vector<Vcontrolled_Isource*>& G,
	std::vector<iSource_nonlin*>& Isource_non_lin, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv, const std::vector<Voltage_source*> Vsourcesdc,
	const std::vector<Current_source*> Isourcesdc);
void ac_analysis_ssem(double omega, std::vector<Node*>& nodes, std::vector<ConductanceDevice*>& devices_small_sig, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv,
	std::vector<Vcontrolled_Isource*>& G, const std::vector<Voltage_source*> Vsourcesac, const std::vector<Current_source*>Isourcesac);
#endif



