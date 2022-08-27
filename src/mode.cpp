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

std::string mode_diode(std::complex<double> Va, std::complex<double> Vc) {
	std::string mode;
	double VA = abs(Va);
	double VC = abs(Vc);
	if (VA >= VC) {
		mode = "on";
	}
	else {
		mode = "off";
	}
	return mode;
}
std::string mode_bjt(std::string type, std::complex<double> Vc, std::complex<double> Ve, std::complex<double> Vb) {
	std::string mode;
	double VB = abs(Vb);
	double VC = abs(Vc);
	double VE = abs(Ve);
	double VBC = VB - VC;
	double VBE = VB - VE;
	if (type == "NPN"){
		if (VBC <= 0 && VBE >= 0) {
			mode = "active";
		}
		else if (VBC <= 0 && VBE < 0) {
			mode = "off";
		}
		else if (VBC > 0 && VBE > 0) {
			mode = "saturation";
		}
		else if (VBC > 0 && VBE <= 0) {
			mode = "reverse";
		}
	}
	else if (type == "PNP") {
		if (VBE <= 0 && VBC >= 0) {
			mode = "active";
		}
		else if(VBE<=0, VBC<0){
			mode = "saturation";
		}
		else if (VBE> 0, VBC <= 0) {
			mode = "reverse";
		}
		else {
			mode = "off";
		}
	}
	return mode;
}
std::string mode_mosfet(std::string type, std::complex<double> Vd, std::complex<double> Vg, std::complex<double> Vs, double Vt) {
	std::string mode;
	double VD = abs(Vd);
	double VG = abs(Vg);
	double VS = abs(Vs);
	double VDS = VD - VS;
	double VGS = VG - VS;
	if (type == "NMOS") {
		if (VGS < Vt) {
			mode = "off";
		}
		else if (VGS >= Vt && VDS < VGS - Vt) {
			mode = "linear";
		}
		else if (VGS >= Vt && VDS >= VGS - Vt) {
			mode = "saturation";
		}
	}
	else if (type == "PMOS") {
		if (VGS > Vt) {
			mode = "off";
		}
		else if (VGS <= Vt && VDS > VGS - Vt) {
			mode = "linear";
		}
		else if (VGS <= Vt && VDS <= VGS - Vt) {
			mode = "saturation";
		}
	}
	
	return mode;
}

