#ifndef MONTECARLO_H
#define MONTECARLO_H

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
#include "./ProcessGenerator.h"

typedef double (*Payoff)(ProcessGenerator &process);
typedef std::vector<double> ordinates;

class MonteCarlo{
public:
	MonteCarlo(ProcessGenerator &process, Payoff payoff):
	 process(process), payoff(payoff){}

	ProcessGenerator &process;
	Payoff payoff;

	double expectation(const int nb_trials)
	{
		double cumsum_draw = 0;
		process.generate_paths(nb_trials);

		for (int i = 0; i < nb_trials; i++)
		{
			cumsum_draw += payoff(process);
		}
		return cumsum_draw/nb_trials;
	}

private:
	std::vector<double> draws;

};

#endif
