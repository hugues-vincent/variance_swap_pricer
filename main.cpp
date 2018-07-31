//  runn with:
//  g++ -std=c++11 main.cpp -lboost_iostreams -lboost_system -lboost_filesystem && ./a.out   
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
#include "headers/VarSwapAnalytical.h"


// typedef double (*Payoff)(StochasticModel model_draw);

typedef std::vector<std::pair<double, double>> curve;
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
    nb_trials = 1000;
    // number of steps & Maturity
    N = 1000;
    T = 1;
    // inital conditions for the Heston model SDE
    S_0 = 100;
    V_0 = 50;
    //heston model parameters
    rate = 0.03; // risk free rate
    kappa = 3; // mean reversion coeff
    theta = 150; // long term variance
    sigma = 5; // vol of vol
    rho = 0.5; // correlation
    // TG Scheme parameters
    gamma1 = 0.5; 
    gamma2 = 0.5;

    // display paramters
    std::vector<string> names = std::vector<string>();
    int nb_mc = 2;

    HestonEuler hestonEuler = HestonEuler(T,N, rate, kappa, theta, sigma, rho, S_0, V_0);
    HestonBroadieKayaTG hestonTG = HestonBroadieKayaTG(T,N, rate, kappa, theta, sigma, rho, S_0, V_0, gamma1, gamma2);
    HestonBroadieKayaQE hestonQE = HestonBroadieKayaQE(T,N, rate, kappa, theta, sigma, rho, S_0, V_0, gamma1, gamma2);


    // ===============================
    // Plot model path
    // ===============================
    // std::vector<std::vector<ordinates>> c = { 
    //         hestonEuler.generate_paths(nb_mc, "V"), 
    //         hestonTG.generate_paths(nb_mc, "V"),
    //         hestonQE.generate_paths(nb_mc, "V")
    //     };
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     names.push_back("euler");
    // }
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     names.push_back("tg");
    // }
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     names.push_back("qe");
    // }
    // plot(concat<ordinates>(c), T, names,  hestonEuler.param_to_string());


    // ===============================
    // Curve function of a paramater
    // ===============================
    // curve mc_euler = curve();
    // curve mc_tg = curve();
    // std::vector<double> Ss = {0.1, 0.3, 0.4};
    // for (double S : Ss){
    //     hestonEuler.set_parameters(rate, kappa, theta, S, rho, S_0, V_0);
    //     mc_euler.push_back(make_pair(S, monte_carlo(hestonEuler, nb_trials) * 100));
    //     hestonTG.set_parameters(rate, kappa, theta, S, rho, S_0, V_0);
    //     mc_tg.push_back(make_pair(S, monte_carlo(hestonTG, nb_trials) * 100));

    // }
    // plot({mc_euler, mc_tg}, {"euler", "TG"});

    // ===============================
    // Monte Carlo convergence
    // ===============================
    // vector<curve> mc_paths = vector<curve>();
    // vector<double> path;
    // names = std::vector<string>();
    // curve new_path;
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     new_path = curve();
    //     path = monte_carlo_path(hestonEuler, nb_trials);
    //     for (int i=0 ; i<nb_trials; i++)
    //     {
    //         if( i > 10)
    //         new_path.push_back(make_pair(i, exp(-rate*T) * path[i] * 100));
    //     }
    //     mc_paths.push_back(new_path);
    //     names.push_back("euler");
    // }
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     new_path = curve();
    //     path = monte_carlo_path(hestonTG, nb_trials);
    //     for (int i=0 ; i<nb_trials; i++)
    //     {
    //         if( i > 10)
    //         new_path.push_back(make_pair(i, exp(-rate*T) * path[i] * 100));
    //     }
    //     mc_paths.push_back(new_path);
    //     names.push_back("tg");
    // }
    // for (int j=0 ; j<nb_mc; j++)
    // {
    //     new_path = curve();
    //     path = monte_carlo_path(hestonQE, nb_trials);
    //     for (int i=0 ; i<nb_trials; i++)
    //     {
    //         if( i > 10)
    //         new_path.push_back(make_pair(i, exp(-rate*T) * path[i] * 100));
    //     }
    //     mc_paths.push_back(new_path);
    //     names.push_back("qe");
    // }
    // plot(mc_paths, names, hestonEuler.param_to_string());


    // ===============================
    // print Monte Carlo simple result
    // ===============================  
    
    // double tg = monte_carlo(hestonTG, nb_trials);
    // print("MC hestonTG", tg * 100 );

    // double mc_qe = monte_carlo(hestonQE, nb_trials);
    // print("MC hestonQE", mc_qe * 100);

    VarSwapAnalytical varswap = VarSwapAnalytical(T,N, rate, kappa, theta, sigma, rho, S_0, V_0);
    double analytic =  varswap.var_swap_analytical();
    double euler = monte_carlo(hestonEuler, nb_trials);
    
    print("Analytic",analytic);
    print("MC hestonEuler", euler * 100);

    return 0;
} 
#endif
