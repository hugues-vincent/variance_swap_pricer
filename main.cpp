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
#include "headers/HestonBroadieKayaQE.h"
#include "headers/MonteCarlo.h"


// typedef double (*Payoff)(StochasticModel model_draw);

typedef std::vector<double> ordinates;
int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));
    
    int N, nb_trials; 
    double T; 
    double S_0, V_0; 
	double rate, kappa, theta, sigma, rho; 
    double gamma1, gamma2;
    // number of draw in the monte carlo method
    nb_trials = pow(10,3);
    // number of steps & Maturity
    N = 1000;
    T = 1;
    // inital conditions for the Heston model SDE
    S_0 = 30;
    V_0 = 1;
    //heston model parameters
    rate = 0.05; // risk free rate
    kappa = 0.5; // mean reversion coeff
    theta = 1.4; // long term variance
    sigma = 0.3; // vol of vol
    rho = - 0.6; // correlation
    // TG Scheme parameters
    gamma1 = 0.5; 
    gamma2 = 0.5;


    HestonEuler hestonEuler = HestonEuler(T,N, rate, kappa, theta, sigma, rho, S_0, V_0);
    HestonBroadieKayaTG hestonTG = HestonBroadieKayaTG(T,N, rate, kappa, theta, sigma, rho, S_0, V_0, gamma1, gamma2);
    HestonBroadieKayaQE hestonQE = HestonBroadieKayaQE(T,N, rate, kappa, theta, sigma, rho, S_0, V_0, gamma1, gamma2);

    // plot(hestonEuler.generate_paths(20, "W_v"), T);
    // plot(hestonEuler.generate_paths(20, "V"), T);
    // plot(hestonEuler.generate_paths(20, "S"), T);
    // plot(hestonEuler.generate_paths(20, "LnS"), T);
    // plot(hestonEuler.generate_paths(), T);
    // plot(hestonTG.generate_paths(), T);
    // plot({hestonEuler.V, hestonTG.V}, T, { "hestonEuler", "hestonBroadieKaya"});
    // plot({hestonEuler.lnS, hestonTG.lnS}, T, { "hestonEuler", "hestonBroadieKaya"});
    // plot({hestonQE.lnS}, T);

    std::vector<std::vector<ordinates>> c = { 
            hestonEuler.generate_paths(2), 
            hestonTG.generate_paths(2),
            hestonQE.generate_paths(2)
        };
    plot(concat<ordinates>(c), T, {"euler", "euler", "TG", "TG", "QE", "QE"});


    // double mc_euler = monte_carlo(hestonEuler, nb_trials);
    // double mc_tg = monte_carlo(hestonTG, nb_trials);
    // double mc_qe = monte_carlo(hestonQE, nb_trials);

    // p("MC hestonEuler", mc_euler);
    // p("MC hestonTG", mc_tg);
    // p("MC hestonQE", mc_qe);

    return 0;
} 
#endif
