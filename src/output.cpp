#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include <cmath>
#include"class.h"
#include <complex>

void to_f(std::vector<double>& frequencyv, double fs, double fe, int interval) {
	double f = fs;
	while (f < 1000) {
		frequencyv.push_back(f);
		f = f + 10;
	}
	while (f < fe) {
		f = f + 1000;
		frequencyv.push_back(f);
	}
}
void to_m_p(std::vector<double>& magnitudes, std::vector<double>& phases, const std::vector<std::complex<double>> voltages) {
	double magnitude;
	double phase;
	for (int i = 0; i < voltages.size(); i++) {
		magnitude = abs(voltages[i]);
		phase = arg(voltages[i]);
		magnitudes.push_back(magnitude);
		phases.push_back(phase);
	}
}
