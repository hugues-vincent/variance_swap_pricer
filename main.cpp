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
#include "headers/HestonEulerMaruyama.h"
#include "headers/Payoff.h"
#include "headers/MonteCarlo.h"


// typedef double (*Payoff)(StochasticModel model_draw);


int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));
    
    int N; // number of steps
    double T; // Maturity
    double S_0, V_0; // inital conditions for the Heston model SDE
	double rate, kappa, theta, sigma, rho; //heston model 
    N = 1000;
    T = 1;
    S_0 = 10;
    V_0 = 1;
    rate = 0.05;
    kappa = 10;
    theta = 0.1;
    sigma = 0.3;
    rho = - 0.6;


    HestonEulerMaruyama heston = HestonEulerMaruyama(T,N, rate, kappa, theta, sigma, rho, S_0, V_0);
    std::vector<std::vector<double>> trials;
    plot(heston.generate_paths("V", 20), T);
    // plot(heston.generate_paths("S", 20), T);
    // plot(heston.generate_paths("LnS", 20), T);
    return 0;
} 
#endif
