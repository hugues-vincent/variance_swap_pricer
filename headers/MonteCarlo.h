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

typedef std::vector<double> ordinates;

double monte_carlo(ProcessGenerator &process, const int nb_trials)
{
	double expectation = 0;
	int i25 = int(nb_trials * .25) , i50 = int(nb_trials * .50), i75 = int(nb_trials * .75);
	cout << "==================== " << endl;
	for(int i(0) ; i<nb_trials ; i++)
	{
		if (i == i25  || i == i50 || i == i75)
			cout << "step: " << i << "/" << nb_trials << endl;
		ordinates trial = process.new_trial();
		for(int j(1) ; j < process.N ;  j++)
		{
			expectation += pow(trial[j] - trial[j-1],2);
		}
	}
	expectation = process.dt * expectation / nb_trials;
	cout << endl;
	return expectation;
}

std::vector<double> monte_carlo_path(ProcessGenerator &process, const int nb_trials)
{
	std::vector<double> expectation = vector<double>();
	expectation.push_back(0);
	double current_trial = 0;
	int i25 = int(nb_trials * .25) , i50 = int(nb_trials * .50), i75 = int(nb_trials * .75);
	cout << "==================== " << endl;
	for(int i(0) ; i<nb_trials ; i++)
	{
		current_trial = 0;
		if (i == i25  || i == i50 || i == i75)
			cout << "step: " << i << "/" << nb_trials << endl;
		ordinates trial = process.new_trial();
		for(int j(1) ; j < process.N ;  j++)
		{
			current_trial += pow(trial[j] - trial[j-1],2) ;
		}
		expectation.push_back(process.dt * current_trial / (i+1) + expectation[i] * i / (i+1));
	}
	cout << endl;
	return expectation;
}
#endif
