#ifndef STOCHASTICMODEL_H
#define STOCHASTICMODEL_H

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

#include "utils/plot.h"
#include "utils/utils.h"

class StochasticModel {
public:
	StochasticModel(double T, int N):T(T), N(N), dt(T/N), S(N){}
    std::vector<double> S;
    double T, N, dt;
    virtual void new_trial()
    {
    	S = std::vector<double>(N);
    }
};

class Heston: public StochasticModel {
public:
    Heston(double T, int N, double r, double k, double t, double s, double rho):
    StochasticModel(T,N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), V(N){
    	new_trial();
    }
    std::vector<double> W_s, W_v, V;
    double rate, kappa, theta, sigma, rho;

 	void new_trial() override
    {
    	correlated_draws(W_v, W_s);
    	calc_vol_path(W_v, V);
    	calc_spot_path(W_s, V, S);
    }
private:
	void correlated_draws(vector<double>& uncorr_draws, vector<double>& corr_draws)
	{
		std::vector<double> delta_uncorr_draws(N), delta_corr_draws(N);
	    for (int i(1); i<corr_draws.size(); i++)
	    {
			delta_uncorr_draws[i] = sqrt(dt) * gaussian_draw();
			delta_corr_draws[i] = rho * delta_uncorr_draws[i] + sqrt(1 - rho*rho) * sqrt(dt) * gaussian_draw();
			uncorr_draws[i] = uncorr_draws[i-1] + delta_uncorr_draws[i];
			corr_draws[i] = corr_draws[i-1] + delta_corr_draws[i];
	    }
	}
	void calc_vol_path(const vector<double>& vol_draws, vector<double>& vol_path)
	{
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        vol_path[i] = vol_path[i-1] + kappa * (theta - v_max) * dt + sigma * sqrt(v_max * dt) * vol_draws[i-1];
	    }
	}

	void calc_spot_path(const vector<double>& spot_draws, const vector<double>& vol_path, vector<double>& spot_path)
	{	    
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        spot_path[i] = spot_path[i-1] * exp( (rate - 0.5*v_max))*dt + sqrt(v_max*dt)*spot_draws[i-1];
	    }
	}

};

#endif
