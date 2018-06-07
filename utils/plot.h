#ifndef PLOT_H
#define PLOT_H
#include "../plugins/gnuplot-iostream.h"
#include <string>
#include "../headers/StochasticModel.h"

using namespace std;
using curve = std::vector<std::pair<double, double>>;

void plot(vector<pair<double, double>> xy_pts, const double x_l, const double x_r, const double y_l, const double y_r, const string name="")  
{
    	Gnuplot gp;

		vector<pair<double, double> > xy_pts_A;
		for(double x=-2; x<2; x+=0.01) {
			double y = x*x*x;
			xy_pts_A.push_back(make_pair(x, y));
		}

		vector<pair<double, double> > xy_pts_B;
		for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
			double theta = alpha*2.0*3.14159;
			xy_pts_B.push_back(make_pair(cos(theta), sin(theta)));
		}

		gp << "set xrange ["<<x_l<<":"<< x_r<<"]\nset yrange ["<<y_l<<":"<< y_r<<"]\n" << "plot";
		gp << gp.file1d(xy_pts) <<  "with lines title '"<< name <<"',"<< endl;

	#ifdef _WIN32
		// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
		// the gnuplot window doesn't get closed.
		cout << "Press enter to exit." << endl;
		cin.get();
	#endif 
};
void plot(std::vector<curve> xy_ptss, const double x_l, const double x_r, const double y_l, const double y_r, vector<string> names = vector<string>()) 
{
    	Gnuplot gp;
    	int N = xy_ptss.size();
    	if(names.empty())
    		names =  vector<string>(N, "");
		gp << "set xrange ["<<x_l<<":"<< x_r<<"]\n set yrange ["<<y_l<<":"<< y_r<<"]\n" << "plot";
		for(int i(0); i<N; i++)
		{
			gp << gp.file1d(xy_ptss[i]) <<  "with lines title '"<< names[i]<<"',";			
		}
		gp << endl;			

};

void plot_stochastic_model(StochasticModel &stochasticmodel){

	int N = stochasticmodel.N;
	double T = stochasticmodel.T;
	double max_= 0 , min_= 0;
	curve model_curve;
    std::vector<double> xs = std::vector<double>(N);

    
    for(int i(1); i<N; i++)
    {
    	xs[i] = xs[i-1] + T/N;
    }
	for(double i=0; i<N; i++) 
	{
		model_curve.push_back(std::make_pair(xs[i], stochasticmodel.S[i]));
		max_ = max(stochasticmodel.S[i],max_);
		min_ = min(stochasticmodel.S[i],min_);
	}
	max_ += 0.1;
	min_ -= 0.1;
    plot(model_curve, 0, T, min_, max_);
}
void plot_stochastic_model(StochasticModel &stochasticmodel, int nb_trials){

	int N = stochasticmodel.N;
	double max_= 0 , min_= 0, T = stochasticmodel.T;
	std::vector<curve> model_curves;
    std::vector<double> xs = std::vector<double>(N);
	curve model_curve;
    
    for(int i(1); i<N; i++)
    	xs[i] = xs[i-1] + T/N;
	
	for(int i(0); i<nb_trials; i++)
	{
		stochasticmodel.new_trial();
		model_curve = curve();
		for(double i=0; i<N; i++) {
			model_curve.push_back(std::make_pair(xs[i], stochasticmodel.S[i]));
			max_ = max(stochasticmodel.S[i],max_);
			min_ = min(stochasticmodel.S[i],min_);
		}
		model_curves.push_back(model_curve);
	}
	max_ += 0.1;
	min_ -= 0.1;
    plot(model_curves, 0, T, min_, max_);
}
#endif
