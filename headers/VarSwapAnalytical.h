#ifndef VAR_SWAP_ANALYTICAL_H
#define VAR_SWAP_ANALYTICAL_H

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
#include <complex>

#include "../utils/plot.h"
#include "../utils/utils.h"
#include "./ProcessGenerator.h"
#include "./HestonModel.h"

typedef std::vector<double> ordinates;
typedef complex<double> cd;

double C(double w, double tau, double rate, double kappa, double theta, double sigma, double rho)
{
	double a = kappa - rho * sigma * w;
	double b =0 ;
	double g = (a + b) / (a - b);
	return rate * (w - 1) * tau + kappa*theta/sigma2 * ((a + b) * tau - 2 * ln((1 - g*exp(b*tau)) / 1 - g)); 
}

double var_swap_analytical(double T, int N, double rate, double kappa, double theta, double sigma, double rho, double S_0, double V_0) 
{
	double tau, w;
	double c_i, w_i, t_i;
	double a, b, g, q;
	double C1, C2, D1, D2;
	double g_v0, sum_gi_vo;
	double sigma2 = pow(sigma, 2);
	
	tau =  T/N;
	w = 0;
	q = 2 * kappa * theta / sigma2;
	a = kappa - rho * sigma * w;
	b = sqrt(pow(a, 2) + sigma2 * (w*w  + w));
	g = (a + b) / (a - b);

	C = rate * (w - 1) * tau + kappa*theta/sigma2 * ((a + b) * tau - 2 * ln((1 - g*exp(b*tau)) / 1 - g)); 
	D = (a + b) * (1 - exp(b*tau)) / sigma2 / (1 - g*exp(b*tau)); 


	g_v0 = 0 ;

	sum_gi_vo = 0;
	for (int i(2) ; i<N ; i++) 
	{
		t_i = (i -1) * tau ;
		c_i = 2 * kappa / sigma2 / (1 - exp(-kappa*t_i));
		w_i = c_i * V_0 * exp(-kappa*t_i);
		sum_gi_vo += pow(D1, 2) * (q + 2*w_i + pow(q + w_i, 2))/pow(c_i, 2) + (2*C1*D1 - D2)*(q + w_i)/c_i + pow(C1, 2) - C2;
	}
	
	return pow(100, 2)(g_vo + sum_gi_vo)/T;
}
#endif
