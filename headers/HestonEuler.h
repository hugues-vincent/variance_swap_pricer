#ifndef HESTON_EULER_H
#define HESTON_EULER_H

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

class HestonEuler : public ProcessGenerator, public HestonModel {
public:
    HestonEuler(double T, int N, double r, double k, double t, double s, double rho, double S_0, double V_0):
    ProcessGenerator(N, T), HestonModel(r, k, t, s, rho, S_0, V_0){
    	new_trial();
    }

    ordinates W_s, W_v;

 	ordinates new_trial()
    {
		reset_paths();
    	correlated_draws(W_v, W_s);
    	generate_vol_path(W_v, V);
    	generate_spot_path(W_s, V, S);
    	generate_log_spot_path(W_s, V, lnS);
    	return lnS;
    }
    std::vector<ordinates> generate_paths(const int nb_trials = 1, const string path_name = "lnS")
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
	void correlated_draws(ordinates& uncorr_draws, ordinates& corr_draws)
	{
		ordinates delta_uncorr_draws(N), delta_corr_draws(N);
	    for (int i(1); i<corr_draws.size(); i++)
	    {
			delta_uncorr_draws[i] = sqrt(dt) * gaussian_draw();
			delta_corr_draws[i] = rho * delta_uncorr_draws[i] + sqrt(1 - rho*rho) * sqrt(dt) * gaussian_draw();
			uncorr_draws[i] = uncorr_draws[i-1] + delta_uncorr_draws[i];
			corr_draws[i] = corr_draws[i-1] + delta_corr_draws[i];
	    }
	}
	void generate_vol_path(const ordinates& vol_draws, ordinates& vol_path)
	{
		vol_path[0] = V_0;
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_v = vol_draws[i] - vol_draws[i-1];
	        vol_path[i] = vol_path[i-1] + kappa * (theta - v_max) * dt + sigma * sqrt(v_max) * dW_v;
	    }
	}

	void generate_spot_path(const ordinates& spot_draws, const ordinates& vol_path, ordinates& spot_path)
	{	
		spot_path[0] = S_0;

	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_s = spot_draws[i] - spot_draws[i-1];
	        spot_path[i] = spot_path[i-1] + rate*spot_path[i-1]*dt + sqrt(v_max)*spot_path[i-1]*dW_s;
	    }
	}
	void generate_log_spot_path(const ordinates& spot_draws, const ordinates& vol_path, ordinates& log_spot_path)
	{	
		
		log_spot_path[0] = log(S_0);
	    for (int i(1); i<N; i++)
	    {
	        double v_max = max(vol_path[i-1], 0.0);
	        double dW_s = spot_draws[i] - spot_draws[i-1];
	        log_spot_path[i] = log_spot_path[i-1] + (rate - v_max/2)*dt + sqrt(v_max)*dW_s;
	    }
	}

	void reset_paths()
	{
		V.clear();
		S.clear();
		lnS.clear();
		W_s.clear();
		W_v.clear();
		for (int i(0); i<N; i++)
		{
			V.push_back(0);
			S.push_back(0);
			lnS.push_back(0);
			W_s.push_back(0);
			W_v.push_back(0);
		}
	}
};

#endif
