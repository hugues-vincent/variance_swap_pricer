#ifndef HESTON_BROADIE_KAYA_QE_H
#define HESTON_BROADIE_KAYA_QE_H

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
#include "./HestonModel.h"

typedef std::vector<double> ordinates;

class HestonBroadieKayaQE: public ProcessGenerator, public HestonModel {
public:
    HestonBroadieKayaQE(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0, double gamma1, double gamma2): 
    ProcessGenerator(N, T), HestonModel(r, k, t, s, rho, S_0, V_0), gamma1(gamma1), gamma2(gamma2){
    	k0 = -rho * kappa * theta / sigma * dt;
    	k1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho /sigma;
    	k2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho /sigma;
    	k3 = gamma1 * (1 - pow(rho,2)) * dt;
    	k4 = gamma2 * (1 - pow(rho,2)) * dt;
		exp_kapa = exp(-kappa * dt);
    	new_trial();
    }
    double gamma1, gamma2;
    double exp_kapa, alpha;

 	ordinates new_trial()
    {
		reset_paths();
		generate_vol_path(V);
		generate_log_spot_path(V, lnS);
		return lnS;
    }
	void set_parameters(double r_, double k_, double t_, double s_, double rho_, double S_0_, double V_0_)
	{
    	rate = r_;
    	kappa = k_;
    	theta = t_;
    	sigma = s_;
    	rho = rho_;
    	S_0 = S_0_;
    	V_0 = V_0_;
    	k0 = -rho * kappa * theta / sigma * dt;
    	k1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho /sigma;
    	k2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho /sigma;
    	k3 = gamma1 * (1 - pow(rho,2)) * dt;
    	k4 = gamma2 * (1 - pow(rho,2)) * dt;
		exp_kapa = exp(-kappa * dt);
    	new_trial();
	}
    std::vector<ordinates> generate_paths(const int nb_trials = 1, const string path_name = "lnS")
    {
	    std::vector<ordinates> trials;
		for(int i(0); i<nb_trials; i++)
		{
		    new_trial();
		    if (path_name == "V")
			    trials.push_back(V);
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
		double m, s2, psi; 
		double a, b, p, beta; 
		double u;
		int print_limit = 4;
		vol_path[0] = V_0;

	    for (int i(1); i<N; i++)
	    {
	    	m = theta + (vol_path[i-1] - theta) * exp_kapa;
	    	s2 = vol_path[i-1] * pow(sigma, 2) * exp_kapa * (1 - exp_kapa) / kappa  + theta * pow(sigma, 2) * pow((1 - exp_kapa), 2) / (2 * kappa);
	    	
	    	psi =s2/pow(m, 2) ;
			if (0 < psi <= 1.5 ) 
			{
				b = sqrt(2/psi - 1 + sqrt(2/psi) * sqrt(2/psi - 1));
				a = m /(1 + b*b);
		        vol_path[i] =  a * pow(b + gaussian_draw(), 2);
			}
			else 
			{
				p = (psi - 1) / (psi + 1);
				beta = 2 / (m*(psi + 1));
		        u = uniform_draw();
		        if (u > p)
		        	vol_path[i] = log((1 - p) / (1 - u)) / beta;
		        else 
		        	vol_path[i] = 0;
			}
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
	void reset_paths()
	{
		V.clear();
		S.clear();
		lnS.clear();
		for (int i(0); i<N; i++)
		{
			V.push_back(0);
			S.push_back(0);
			lnS.push_back(0);
		}
	}
};

#endif
