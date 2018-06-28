#ifndef PLOT_H
#define PLOT_H
#include "../plugins/gnuplot-iostream.h"
#include <string>
#include "../headers/StochasticModel.h"

using namespace std;
using curve = std::vector<std::pair<double, double>>;

curve path_to_curve(const std::vector<double> path, const double T) 
{
	curve model_curve = curve();
	int N = path.size();
	double x = 0, dt = T/N;
	for(double j=0; j<N; j++) {
		model_curve.push_back(std::make_pair(x, path[j]));
		x += dt;
	}
	return model_curve;
}

void plot(std::vector<std::vector<double>> paths, const double T, vector<string> names = vector<string>()) 
{
	double max_= 0 , min_= 0;
	Gnuplot gp;
	int nb_curve = paths.size();
	for(int i(0); i<nb_curve; i++)
	{
		for(double j=0; j<paths[i].size(); j++) {
			max_ = max(paths[i][j],max_);
			min_ = min(paths[i][j],min_);
		}
	}
	max_ += 0.1;
	min_ -= 0.1;

	if(names.empty())
		names =  vector<string>(nb_curve, "");
	gp << "set xrange ["<<0<<":"<< T<<"]\n set yrange ["<<min_<<":"<< max_<<"]\n" << "plot";
	for(int i(0); i<nb_curve; i++)
	{
		gp << gp.file1d(path_to_curve(paths[i], T)) <<  "with lines title '"<< names[i]<<"',";	
	}
	gp << endl;			

};
#endif
