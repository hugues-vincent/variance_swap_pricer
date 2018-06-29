#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cmath>
#include <cstdlib> //rand
#include <time.h>
#include <fstream>
#include <string>
#include <random>

#include <numeric>
#include <vector>
#include <iostream>
#include <iterator>
#include <functional>

using namespace std;
std::random_device rd{};
std::mt19937 gen{rd()};

double max(const double x, const double y){ return x>y? x : y;}
double min(const double x, const double y){ return x<y? x : y;}

double uniform_draw(const double a, const double b)
{
    std::uniform_int_distribution<> dis{0, 99};
	return  dis(gen);
}

double gaussian_draw(const double  mu = 0, const double  sigma = 1)
{
    std::normal_distribution<> dis{mu, sigma};
    return dis(gen);}

double gaussian_cdf_inverse(const double  u)
{
	// Returns the inverse of cumulative normal distribution function.
	// Reference: Moro, B., 1995, "The Full Monte," RISK (February), 57-58.
	double a[4] = {
		2.50662823884,
		-18.61500062529,
		41.39119773534,
		-25.44106049637
	};

	double b[4] = {
		-8.47351093090,
		23.08336743743,
		-21.06224101826,
		3.13082909833
	};

	double c[9] = {
		0.3374754822726147,
		0.9761690190917186,
		0.1607979714918209,
		0.0276438810333863,
		0.0038405729373609,
		0.0003951896511919,
		0.0000321767881768,
		0.0000002888167364,
		0.0000003960315187
	};

	double x, r;

	x = u - 0.5;
	if( fabs(x) < 0.42 )
	{ 
		r = x * x;
		r = x * ((( a[3]*r + a[2]) * r + a[1]) * r + a[0])/
		((((b[3] * r+ b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
		return (r);
	}

	r = u;
	if( x > 0.0 ) r = 1.0 - u;
	r = log(-log(r));
	r = c[0] + r * (c[1] + r * 
		(c[2] + r * (c[3] + r * 
			(c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r*c[8])))))));
	if( x < 0.0 ) r = -r;

	return (r);
}

void p(const std::vector<double> v)
{
	for(int i(0); i<v.size(); i++)
	{
		cout << v[i] <<" ";
	}
	cout << "\n";
}
void p(const string var_name)
{
	cout << var_name << "\n";
}
void p(const string txt, const double var_name)
{
	cout << txt << ": " << var_name << "\n";
}

void p(const double var_name)
{
	cout << var_name << "\n";
}
void p(const string txt, const int var_name)
{
	cout << txt << ": " << var_name << "\n";
}

void p(const int var_name)
{
	cout << var_name << "\n";
}
#endif
