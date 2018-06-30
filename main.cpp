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
#include "headers/HestonEuler.h"
#include "headers/HestonBroadieKayaTG.h"
#include "headers/MonteCarlo.h"


// typedef double (*Payoff)(StochasticModel model_draw);


int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));
    
    int N; 
    double T; 
    double S_0, V_0; 
	double rate, kappa, theta, sigma, rho; 
    double gamma1, gamma2;
    // number of steps & Maturity
    N = 1000;
    T = 1;
    // inital conditions for the Heston model SDE
    S_0 = 0.04;
    V_0 = 1;
    //heston model parameters
    rate = 0.05;
    kappa = 0.5;
    theta = 0.04;
    sigma = 1;
    rho = - 0.6;
    // TG Scheme parameters
    gamma1 = 0.5; 
    gamma2 = 0.5;


    HestonEuler hestonEuler = HestonEuler(T,N, rate, kappa, theta, sigma, rho, S_0, V_0);
    HestonBroadieKayaTG hestonTG = HestonBroadieKayaTG(T,N, rate, kappa, theta, sigma, rho, S_0, V_0, gamma1, gamma2);

    // plot(hestonEuler.generate_paths("W_v", 20), T);
    // plot(hestonEuler.generate_paths("V", 40), T);
    // plot(heston.generate_paths("S", 20), T);
    // plot(heston.generate_paths("LnS", 20), T);
    

    
    // plot(hestonBroadieKaya.generate_paths("V", 1), T);
    // plot({hestonEuler.V, hestonTG.V}, T, { "hestonEuler", "hestonBroadieKaya"});
    plot({hestonEuler.lnS, hestonTG.lnS}, T, { "hestonEuler", "hestonBroadieKaya"});

    return 0;
} 
#endif
