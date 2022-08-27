#pragma once
#ifndef OUTPUT
#define OUTPUT
#include "class.h"

void to_f(std::vector<double>& frequencyv, double fs, double fe, int interval);
void to_m_p(std::vector<double>& magnitudes, std::vector<double>& phases, const std::vector<std::complex<double>> voltages);

#endif