#ifndef HESTON_BROADIE_KAYA_TG_H
#define HESTON_BROADIE_KAYA_TG_H

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

class HestonBroadieKayaTG: public ProcessGenerator {
public:
    HestonBroadieKayaTG(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0, double gamma1, double gamma2):
    T(T), N(N), dt(T/N), S(N), lnS(N), V(N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), S_0(S_0), V_0(V_0), gamma1(gamma1), gamma2(gamma2){
    	k0 = -rho * kappa * theta / sigma * dt;
    	k1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho /sigma;
    	k2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho /sigma;
    	k3 = gamma1 * (1 - pow(rho,2)) * dt;
    	k4 = gamma2 * (1 - pow(rho,2)) * dt;
    	new_trial();
    }

    ordinates W_s, W_v, V, S, lnS;
    double rate, kappa, theta, sigma, rho;
    double T, N, dt;
    double  S_0, V_0;
    double  gamma1, gamma2;

 	void new_trial()
    {
		generate_vol_path(V);
		generate_log_spot_path(V, lnS);
    }
	std::vector<ordinates> generate_paths(const int nb_trials = 1)
	{
	    std::vector<ordinates> trials;
		for(int i(0); i<nb_trials; i++)
		{
		    new_trial();
				trials.push_back(lnS);
		}
		return trials;	
	}
    std::vector<ordinates> generate_paths(const string path_name, const int nb_trials)
    {
	    std::vector<ordinates> trials;
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
	double k0, k1, k2, k3, k4;
	
	void generate_vol_path(ordinates& vol_path)
	{
		double m, s2, phi; 
		double exp_kapa = exp(-kappa * dt);

		vol_path[0] = V_0;
	    for (int i(1); i<N; i++)
	    {
	    	m = theta + (vol_path[i-1] - theta) * exp_kapa;
	    	s2 = vol_path[i-1] * pow(sigma, 2) * exp_kapa / kappa * (1 - exp_kapa) + theta * pow(sigma, 2)  / (2 * kappa) * pow((1 - exp_kapa), 2);

	        vol_path[i] = max(m + sqrt(s2) * gaussian_draw(), 0);
	    }
	}
	void generate_log_spot_path(const ordinates& vol_path, ordinates& log_spot_path)
	{	
		
		log_spot_path[0] = log(S_0);
	    for (int i(1); i<N; i++)
	    {
	        log_spot_path[i] = log_spot_path[i-1] + k0 + k1 * vol_path[i-1] + k2 * vol_path[i] + sqrt(k3 * vol_path[i-1] + k4 * vol_path[i]) * gaussian_draw();
	    }
	}
};

#endif
