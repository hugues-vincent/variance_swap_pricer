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
void plot(std::vector<std::vector<double>> ys, const double T, vector<string> names = vector<string>()) 
{
	double max_= 0 , min_= 0;
	Gnuplot gp;
	curve model_curve;
	int nb_curve = ys.size();
	for(int i(0); i<nb_curve; i++)
	{
		for(double j=0; j<ys[i].size(); j++) {
			max_ = max(ys[i][j],max_);
			min_ = min(ys[i][j],min_);
		}
	}
	max_ += 0.1;
	min_ -= 0.1;

	if(names.empty())
		names =  vector<string>(nb_curve, "");
	gp << "set xrange ["<<0<<":"<< T<<"]\n set yrange ["<<min_<<":"<< max_<<"]\n" << "plot";
	for(int i(0); i<nb_curve; i++)
	{
		model_curve = curve();
		int N = ys[i].size();
		double x = 0;
		for(double j=0; j<N; j++) {
			model_curve.push_back(std::make_pair(x, ys[i][j]));
			x += T/N;
		}
		gp << gp.file1d(model_curve) <<  "with lines title '"<< names[i]<<"',";	
	}
	gp << endl;			

};
#endif
