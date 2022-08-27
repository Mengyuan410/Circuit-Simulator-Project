#ifndef TO_EQUIV
#define ITERSTION
#include "class.h"

void dc_convert(std::vector <ConductanceDevice*>& large_sig_devices, std::vector<Vcontrolled_Isource*>& G, std::vector<iSource_nonlin*>& isource_non_lin, std::vector<Resistor_non_lin_equiv*>& non_lin_resis_equiv, std::vector<ConductanceDevice*> resistors, std::vector<BJT*> BJTs, std::vector<Mosfet*> Mosfets,
    std::vector<Diode*> Diodes);
void ac_convert(std::vector<ConductanceDevice*>& ssem_devices, const std::vector<Resistor_non_lin_equiv*> non_lin_resis_equiv, const std::vector<ConductanceDevice*> rlcdevices);
#endif




