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
#include "headers/MonteCarlo.h"


typedef double (*Payoff)(StochasticModel model_draw);


int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));

	std::vector<std::pair<double, double> > xy_pts_B;
	for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
		double theta = alpha*2.0*3.14159;
		xy_pts_B.push_back(std::make_pair(cos(theta), sin(theta)));
	}
    
    int N(10000); // Number of steps
    double T(1.0); // Maturity
    double S_0(100.0), v_0(0.5); // inital conditions for the Heston model SDE
	double rate(0.5), kappa(0.3), theta(.9), sigma(0.9), rho(0.6), tau(1.0/8.0); //heston model parameters

    Heston hm = Heston(T,N, rate, kappa, theta, sigma, rho);
    plot_stochastic_model(hm, 7);
    return 0;
} 
#endif
