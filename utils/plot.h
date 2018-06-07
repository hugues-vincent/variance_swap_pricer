#include "../plugins/gnuplot-iostream.h"
#include <string>

using namespace std;

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
void plot(vector<vector<pair<double, double>>> xy_ptss, const double x_l, const double x_r, const double y_l, const double y_r, vector<string> names = vector<string>()) 
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


	// #ifdef _WIN32
	// 	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// 	// the gnuplot window doesn't get closed.
	// 	cout << "Press enter to exit." << endl;
	// 	cin.get();
	// #endif 
};
