//  runn with:
//  g++ -std=c++11 hv.cpp -lboost_iostreams -lboost_system -lboost_filesystem && ./a.out   
#ifndef MAIN_CPP
#define MAIN_CPP

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
#include "headers/StochasticModel.h"
#include "headers/Payoff.h"


typedef double (*Payoff)(StochasticModel model_draw);

class MonteCarlo{
public:
	MonteCarlo(int n, StochasticModel model_draw, Payoff payoff):
	 model_draw(model_draw), nb_trials(n), payoff(payoff){}
	int nb_trials;
	StochasticModel model_draw;
	Payoff payoff;

	double expectation()
	{
		double cumsum_draw = 0;
		for (int i = 0; i < nb_trials; i++)
		{
			model_draw.new_trial();
			cumsum_draw += payoff(model_draw);
		}
		return cumsum_draw/nb_trials;
	}

private:
	std::vector<double> draws;

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
    
    double rate(0.5) ,kappa(0.3) ,theta(.9) ,sigma(0.9) ,rho(0.6), tau(1.0/8.0);

    Heston hm = Heston(T,N, rate, kappa, theta, sigma, rho);


    // plot_stochastic_model(hm, 7);
    return 0;
} 
#endif
