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


// typedef double (*Payoff)(StochasticModel model_draw);


int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));
    
    int N(1000); // Number of steps
    double T(1.0); // Maturity
    double S_0(100.0), v_0(0.5); // inital conditions for the Heston model SDE
	double rate(0.5), kappa(0.3), theta(.9), sigma(0.9), rho(0.6), tau(1.0/8.0); //heston model parameters

    Heston heston = Heston(T,N, rate, kappa, theta, sigma, rho);
    // MonteCarlo(100, heston, (Payoff)varswap)
    // plot_stochastic_model(heston, 7);
    plot({heston.W_s, heston.W_v}, T);
    return 0;
} 
#endif
