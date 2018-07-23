#ifndef PLOT_H
#define PLOT_H
#include "../plugins/gnuplot-iostream.h"
#include <string>
#include <float.h>
#include <math.h>

using namespace std;
typedef std::vector<std::pair<double, double>> curve;
typedef std::vector<double> ordinates;

curve path_to_curve(const ordinates path, const double T) 
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

void plot(std::vector<ordinates> paths, const double T, vector<string> names = vector<string>()) 
{
	double max_= -DBL_MAX , min_= DBL_MAX;
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
void plot(std::vector<curve> curve, vector<string> names = vector<string>(), string title = "") 
{
	double max_y= -DBL_MAX , min_y= DBL_MAX, max_x= -DBL_MAX , min_x= DBL_MAX;
	Gnuplot gp;
	int nb_curve = curve.size();
	for(int i(0); i<nb_curve; i++)
	{
		for(double j=0; j<curve[i].size(); j++) {
			max_x = max(curve[i][j].first, max_x);
			min_x = min(curve[i][j].first, min_x);
			max_y = max(curve[i][j].second, max_y);
			min_y = min(curve[i][j].second, min_y);

		}
	}

	if(names.empty())
		names =  vector<string>(nb_curve, "");
	gp << "set title \"" << title << "\" \n set xrange ["<<min_x<<":"<< max_x<<"]\n set yrange ["<<min_y<<":"<< max_y<<"]\n" << " plot ";
	for(int i(0); i<nb_curve; i++)
	{
		gp << gp.file1d(curve[i]) <<  "with lines title '"<< names[i]<<"',";	
	}
	gp << endl;			

};
#endif
