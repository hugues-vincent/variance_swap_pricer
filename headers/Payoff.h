#ifndef PAYOFF_H
#define PAYOFF_H

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
#include "StochasticModel.h"


// double varswap(StochasticModel model_draw)
// {
// 	double cum_log_return = 0;

// 	std::vector<double> S = model_draw.S;
// 	double N = model_draw.N;
// 	double T = model_draw.T;

// 	for(int i(1); i<N; i++)
// 	{
// 		cum_log_return += log(pow(S[i]/S[i-1],2));
// 	}
// 	return pow(100,2)*cum_log_return/T;
// }
#endif
