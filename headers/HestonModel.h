#ifndef HESTON_MODEL_H
#define HESTON_MODEL_H

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

typedef std::vector<double> ordinates;

class HestonModel {
public:
    HestonModel(double r, double k, double t, double s, double rho, double S_0, double V_0): rate(r), kappa(k), theta(t), sigma(s), rho(rho), S_0(S_0), V_0(V_0){}
    double rate, kappa, theta, sigma, rho;
    double  S_0, V_0;
    ordinates V, S, lnS;
};

#endif
