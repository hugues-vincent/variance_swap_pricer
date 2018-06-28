#ifndef HESTON_EULER_MARUYAMA_H
#define HESTON_EULER_MARUYAMA_H

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

class HestonEulerMaruyama {
public:
    HestonEulerMaruyama(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0):
    T(T), N(N), dt(T/N), S(N), LnS(N), V(N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), S_0(S_0), V_0(V_0){
    	new_trial();
    }
    std::vector<double> W_s, W_v, V, S, LnS;
    double rate, kappa, theta, sigma, rho, T, N, dt, S_0, V_0;

 	void new_trial()
    {
    	correlated_draws(W_v, W_s);
    	calc_vol_path(W_v, V);
    	calc_spot_path(W_s, V, S);
    	calc_log_spot_path(W_s, V, LnS);
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
				trials.push_back(LnS);
		}
		return trials;
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
		vol_path[0] = V_0;
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_v = vol_draws[i] - vol_draws[i-1];
	        vol_path[i] = vol_path[i-1] + kappa * (theta - v_max) * dt + sigma * sqrt(v_max) * dW_v;
	    }
	}

	void calc_spot_path(const vector<double>& spot_draws, const vector<double>& vol_path, vector<double>& spot_path)
	{	
		spot_path[0] = S_0;

	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_s = spot_draws[i] - spot_draws[i-1];
	        spot_path[i] = spot_path[i-1] + rate*spot_path[i-1]*dt + sqrt(v_max)*spot_path[i-1]*dW_s;
	    }
	}
	void calc_log_spot_path(const vector<double>& spot_draws, const vector<double>& vol_path, vector<double>& log_spot_path)
	{	
		
		log_spot_path[0] = log(S_0);
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_s = spot_draws[i] - spot_draws[i-1];
	        log_spot_path[i] = log_spot_path[i-1] + (rate - v_max/2)*dt + sqrt(v_max)*dW_s;
	    }
	}

};

#endif
