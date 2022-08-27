#ifndef MODE
#define MODE
#include <string>

std::string mode_bjt(std::string type, std::complex<double> Vc, std::complex<double> Vb, std::complex<double> Ve);
std::string mode_diode(std::complex<double> Va, std::complex<double> Vc);
std::string mode_mosfet(std::string type, std::complex<double> Vd, std::complex<double> Vg, std::complex<double> Vs, double Vt);
#endif
