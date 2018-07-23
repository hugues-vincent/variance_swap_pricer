#ifndef HESTON_MODEL_H
#define HESTON_MODEL_H

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

class HestonModel {
public:
    HestonModel(double r, double k, double t, double s, double rho, double S_0, double V_0): rate(r), kappa(k), theta(t), sigma(s), rho(rho), S_0(S_0), V_0(V_0){}
    double rate, kappa, theta, sigma, rho;
    double  S_0, V_0;
    ordinates V, S, lnS;
    virtual void set_parameters(double r_, double k_, double t_, double s_, double rho_, double S_0_, double V_0_)
    {
    	rate = r_;
    	kappa = k_;
    	theta = t_;
    	sigma = s_;
    	rho = rho_;
    	S_0 = S_0_;
    	V_0 = V_0_;
    }
    string param_to_string()
    {
    	 
    	std::ostringstream oss;
    	oss << "   rate: " << rate 
	    	<< "   kappa: " << kappa 
	    	<< "   theta: " << theta 
	    	<< "   sigma: " << sigma 
	    	<< "   rho: " << rho 
	    	<< "   S_0: " << S_0 
	    	<< "   V_0: " << V_0;
		return oss.str();
    }
};

#endif
