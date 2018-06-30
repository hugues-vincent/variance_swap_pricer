#ifndef PROCESS_GENERATOR_H
#define PROCESS_GENERATOR_H

#include <iostream>
#include <cmath>
#include <cstdlib> //rand
#include <time.h>
#include <fstream>
#include <string>

#include <numeric>
#include <vector>
#include <iostream>
#include <iterator>
#include <functional>

#include "../utils/plot.h"
#include "../utils/utils.h"

typedef std::vector<double> ordinates;

class ProcessGenerator {
public:
	ProcessGenerator(double N, double T): T(T), N(N), dt(T/N) {}
	double N, T, dt;
	virtual std::vector<ordinates> generate_paths(const int nb_trials) = 0;
};

#endif
