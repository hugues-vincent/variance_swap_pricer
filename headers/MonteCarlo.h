#ifndef MONTECARLO_H
#define MONTECARLO_H
#include "StochasticModel.h"

typedef double (*Payoff)(StochasticModel model_draw);
class MonteCarlo{
public:
	MonteCarlo(StochasticModel model_draw, Payoff payoff, int nb_trials):
	 model_draw(model_draw), nb_trials(nb_trials), payoff(payoff){}
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
#endif
