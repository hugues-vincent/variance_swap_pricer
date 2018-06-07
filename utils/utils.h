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

using namespace std;
double max(const double x, const double y){ return x>y? x : y;}
double min(const double x, const double y){ return x<y? x : y;}

double gaussian_draw(double mu = 0, double sigma = 1)
{
    double u = ((double) rand() / (RAND_MAX));
    double v = ((double) rand() / (RAND_MAX));
    return mu + sigma * ( sqrt(-2*log(u)) * sin(2*M_PI*v) );
}
void p(const std::vector<double> v)
{
	for(int i(0); i<v.size(); i++)
	{
		cout << v[i] <<" ";
	}
	cout << "\n";
}
void p(string var_name)
{
	cout << var_name << "\n";
}
void p(double var_name)
{
	cout << std::to_string(var_name) << "\n";
}
void p(int var_name)
{
	cout << std::to_string(var_name) << "\n";
}