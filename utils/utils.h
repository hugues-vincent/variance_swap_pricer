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

// =============================
// Simple utilities function
// =============================
double max(const double x, const double y){ return x>y? x : y;}
double min(const double x, const double y){ return x<y? x : y;}

template <typename T>  
std::vector<T> concat(std::vector<std::vector<T>>& to_concat_vectors)  
{  
	std::vector<T> concat_vector = {};
	for (std::vector<T> &to_concat_vector : to_concat_vectors)
		for (T &to_concat_item : to_concat_vector)
			concat_vector.push_back(to_concat_item);
	return concat_vector;
}  

// =============================
// Random utilities function
// =============================
double uniform_draw(const double a=0, const double b=1)
{
	std::uniform_real_distribution<> dis{a, b};
	return dis(gen);
}

double gaussian_draw(const double  mu = 0, const double  sigma = 1)
{
	std::normal_distribution<> dis{mu, sigma};
	return dis(gen);}


// =============================
// Display utilities function
// =============================
void print(const std::vector<double> v)
{
	for(int i(0); i<v.size(); i++)
	{
		cout << v[i] <<" ";
	}
	cout << "\n";
}
void print(const string var_name)
{
	cout << var_name << "\n";
}
void print(const string txt, const double var_name)
{
	cout << txt << ": " << var_name << "\n";
}

void print(const double var_name)
{
	cout << var_name << "\n";
}
void print(const string txt, const int var_name)
{
	cout << txt << ": " << var_name << "\n";
}

void print(const int var_name)
{
	cout << var_name << "\n";
}
#endif
