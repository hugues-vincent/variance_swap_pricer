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

class VarSwapAnalytical : public HestonModel {
public:
	VarSwapAnalytical(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0) : HestonModel(r, k, t, s, rho, S_0, V_0), T(T), N(N) {}
	double T, N;

	double var_swap_analytical() 
	{
		double tau, q;
		double c_i, w_i, t_i;
		double C1, C2, D1, D2;
		double g_v0, sum_gi_vo;
		double sigma2 = pow(sigma, 2);
		
		tau =  T/N;
		q = 2 * kappa * theta / sigma2;

		C1 = C_derivative_1(tau, pow(10, -5));
		C2 = C_derivative_2(tau, pow(10, -5));
		D1 = D_derivative_1(tau, pow(10, -5));
		D2 = D_derivative_2(tau, pow(10, -5));

		print("C1", C1);
		print("C2", C2);
		print("D1", D1);
		print("D2", D2);

		g_v0 =  pow(D1, 2) * V_0 + (2*C1*D1 - D2)*V_0 + pow(C1, 2) - C2;

		sum_gi_vo = 0;
		for (int i(2) ; i<N ; i++) 
		{
			t_i = i * tau ;
			c_i = 2 * kappa / sigma2 / (1 - exp(-kappa*t_i));
			w_i = c_i * V_0 * exp(-kappa*t_i);
			sum_gi_vo += pow(D1, 2) * (q + 2*w_i + pow(q + w_i, 2))/pow(c_i, 2) + (2*C1*D1 - D2)*(q + w_i)/c_i + pow(C1, 2) - C2;
		}
		// version de l'article
		// return (g_v0 + sum_gi_vo)*100*100/T;
		// version du résumé
		return (g_v0 + sum_gi_vo)*tau;
	}

private:
	double C(double tau, double w)
	{
		double sigma2 = pow(sigma, 2);
		double a = kappa - rho * sigma * w;
		double b = sqrt(pow(a, 2) + sigma2 * (w*w  + w));
		
		// version de l'article
		// double g = (a + b) / (a - b);
		// return rate * (w - 1) * tau + kappa*theta/sigma2 * ((a + b) * tau - 2 * log((1 - g*exp(b*tau)) / (1 - g)));
		
		// version du résumé
		double g = (a - b) / (a + b);
		return rate * (w - 1) * tau + kappa*theta/sigma2 * ((a - b) * tau - 2 * log((1 - g*exp(b*tau)) / (1 - g)));

	}
	double D(double tau, double w)
	{
		double sigma2 = pow(sigma, 2);
		double a = kappa - rho * sigma * w;
		double b = sqrt(pow(a, 2) + sigma2 * (w*w  + w));
		// version de l'article
		// double g = (a + b) / (a - b);
	    // return (a + b) * (1 - exp(b*tau)) / sigma2 / (1 - g*exp(b*tau));

		// version du résumé
		double g = (a - b) / (a + b);
	    return (a - b) * (1 - exp(-b*tau)) / sigma2 / (1 - g*exp(-b*tau));

	}

	double C_derivative_1(double tau, double h)
	{
	    return (C(tau, h) - C(tau, -h))/(2*h);
	}

	double C_derivative_2(double tau, double h)
	{
		print("C(tau, 0)", C(tau, 0));
	    return (C(tau, h) - 2*C(tau, 0) + C(tau, -h)) / (h*h);
	}

	double D_derivative_1(double tau, double h)
	{
	    return (D(tau, h) - D(tau, -h))/(2*h);
	}

	double D_derivative_2(double tau, double h)
	{
	    return (D(tau, h) - 2*D(tau, 0) + D(tau, -h)) / (h*h);
	}
};

#endif
