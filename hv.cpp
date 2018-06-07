//  runn with:
//  g++ -std=c++11 hv.cpp -lboost_iostreams -lboost_system -lboost_filesystem && ./a.out   

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

using curve = std::vector<std::pair<double, double>>;


class PricingModel {
public:
	PricingModel(double T, int N):T(T), N(N), dt(T/N), S(N){}
    std::vector<double> S;
    double T, N, dt;
    virtual void new_trial()
    {
    	S = std::vector<double>(N);
    }
};

typedef void (* payoff)(PricingModel model_draw);

class Heston {
public:
    Heston(double T, int N, double r, double k, double t, double s, double rho):
     T(T), N(N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), V(N), S(N){
    	correlated_draws(W_v, W_s);
    	calc_vol_path(W_v, V);
    	calc_spot_path(W_s, V, S);
     	// new_trial();
    }
    std::vector<double> W_s, W_v, S, V;
    double T, N, rate, kappa, theta, sigma, rho, dt;
     void new_trial() 
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

class MonteCarlo{
public:
	MonteCarlo(int n, double T, int N, PricingModel model_draw):
	 model_draw(T, N),nb_trials(n){}
	
	double expectation()
	{
		double cumsum_draw = 0;
		for (int i = 0; i < nb_trials; i++)
		{
			model_draw.new_trial();

		}
		return 0;
	}
	int nb_trials;
	PricingModel model_draw;

private:
	std::vector<double> draws;
	double varswap(PricingModel model_draw)
	{
		double cum_log_return = 0;

		std::vector<double> S = model_draw.S;
		double N = model_draw.N;
		double T = model_draw.T;

		for(int i(1); i<N; i++)
		{
			cum_log_return += log(pow(S[i]/S[i-1],2));
		}
		return pow(100,2)*cum_log_return/T;
	}
};

int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));

	std::vector<std::pair<double, double> > xy_pts_B;
	for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
		double theta = alpha*2.0*3.14159;
		xy_pts_B.push_back(std::make_pair(cos(theta), sin(theta)));
	}
    
    double T(1.0);
    int N(10000);
    double S_0(100.0);
    double v_0(0.5);
    
    double rate(0.5);
    double kappa(0.3);
    double theta(.9);
    double sigma(0.9);
    double rho(0.6);
    
    double tau(1.0/8.0); // ti - ti-1
    Heston hm = Heston(T,N, rate, kappa, theta, sigma, rho);
    std::vector<double> xs = std::vector<double>(N);
    for(int i(1); i<N; i++)
    	xs[i] = xs[i-1] + T/N;
    p("ws");
    // p(hm.W_s);
    // p("wv");
    // p(hm.W_v);
    // p("V");
    // p(hm.V);
    // p("S");
    // p(hm.S);
    
	std::vector<std::pair<double, double>> ws, wv, s, v, concat;
	int max_= 1 , min_=1;
	for(double i=0; i<N; i++) {
		wv.push_back(std::make_pair(xs[i], hm.W_v[i]));
		ws.push_back(std::make_pair(xs[i], hm.W_s[i]));
		s.push_back(std::make_pair(xs[i], hm.S[i]));
		v.push_back(std::make_pair(xs[i], hm.V[i]));
		max_ = max(hm.W_v[i],max_);
		max_ = max(hm.W_s[i],max_);
		max_ = max(hm.V[i],max_);
		max_ = max(hm.S[i],max_);
		min_ = min(hm.W_v[i],min_);
		min_ = min(hm.W_s[i],min_);
		min_ = min(hm.V[i],min_);
		min_ = min(hm.S[i],min_);

	}
	min_ -= .4;
	max_ += .4;
 	cout << max_ << endl;
 	cout << min_ << endl;
    plot({wv, ws, s, v}, 0, T, -1, 1, {"wv", "ws", "s", "v"});
    // plot(v, 0, T, min_, max_, "v");


    return 0;
} 