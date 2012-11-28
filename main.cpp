#include "lifetime.h"
//#include <fstream>
//#include <gsl/gsl_sf_erf.h>
//#include <cmath>
//#include <cstdlib> //so I can use "exit"
#include <iostream>
//SHOULD ALL THE INCLUDES IN lifetim.cpp BE IN MAIN TOO???

using namespace std;

int main() {

	//const double tau = 0.01;
	//const double sigma = 0.2;
	double meas[10000][2];
	
	read_data(meas);
	//pdf(tau, sigma, meas);
	//nll_tau(meas);
	parabolic_minimiser(meas);
	multimin();
	return 0;
}	
