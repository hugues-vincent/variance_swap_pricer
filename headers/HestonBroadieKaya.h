#ifndef HESTON_BROADIE_KAYA_H
#define HESTON_BROADIE_KAYA_H

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

class HestonBroadieKaya : public ProcessGenerator {
public:
    HestonBroadieKaya(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0):
    T(T), N(N), dt(T/N), S(N), lnS(N), V(N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), S_0(S_0), V_0(V_0){
    	new_trial();
    }
    std::vector<double> W_s, W_v, V, S, lnS;
    double rate, kappa, theta, sigma, rho, T, N, dt, S_0, V_0;

 	void new_trial()
    {
    }
	std::vector<std::vector<double>> generate_paths(const int nb_trials)
	{
	    std::vector<std::vector<double>> trials;
		for(int i(0); i<nb_trials; i++)
		{
		    new_trial();
				trials.push_back(lnS);
		}
		return trials;	
	}
    std::vector<std::vector<double>> generate_paths(const string path_name, const int nb_trials)
    {
	    std::vector<std::vector<double>> trials;
		for(int i(0); i<nb_trials; i++)
		{
		    new_trial();
		    if (path_name == "V")
			    trials.push_back(V);
			else if (path_name == "W_v")
				trials.push_back(W_v);
			else if (path_name == "W_s")
				trials.push_back(W_s);
			else if (path_name == "S")
				trials.push_back(S);
			else
				trials.push_back(lnS);
		}
		return trials;
    }
private:
};

#endif
