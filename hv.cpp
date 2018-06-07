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

#include "./utils/plot.h"
using curve = std::vector<std::pair<double, double>>;

using namespace std;
double gaussian_draw(double mu = 0, double sigma = 1)
{
    double u = ((double) rand() / (RAND_MAX));
    double v = ((double) rand() / (RAND_MAX));
    return mu + sigma * ( sqrt(-2*log(u)) * sin(2*M_PI*v) );
}
void p(const std::vector<double> v)
{
	for(int i(0); i<v.size(); i++)
	{
		cout << v[i] <<" ";
	}
	cout << "\n";
}
void p(string var_name)
{
	cout << var_name << "\n";
}


class Heston {
public:
    Heston(double T, int N, double r, double k, double t, double s, double rho):
    T(T), N(N), rate(r), kappa(k), theta(t), sigma(s), rho(rho), W_s(N), W_v(N), S(N), V(N){
    	correlated_draws(W_v, W_s);
    	calc_vol_path(W_v, V);
    	calc_spot_path(W_s, V, S);

    }
    std::vector<double> W_s, W_v, S, V;
    double T, N, rate, kappa, theta, sigma, rho, dt = T/N;

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
    int N(1000);
    double S_0(100.0);
    double v_0(0.04);
    
    double rate(0.05);
    double kappa(0.3);
    double theta(0.04);
    double sigma(0.9);
    double rho(0.5);
    
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
    
	std::vector<std::pair<double, double> > ws, wv;
	for(double i=0; i<N; i++) {
		wv.push_back(std::make_pair(xs[i], hm.W_v[i]));
		ws.push_back(std::make_pair(xs[i], hm.W_s[i]));
	}
    plot({wv, ws}, 0, 1, -4, 4, {"wv", "ws"});


    return 0;
}