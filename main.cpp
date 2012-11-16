#include "lifetime.h"
//#include <fstream>
//#include <gsl/gsl_sf_erf.h>
//#include <cmath>
//#include <cstdlib> //so I can use "exit"
//#include <iostream>
//SHOULD ALL THE INCLUDES IN lifetim.cpp BE IN MAIN TOO???

using namespace std;

int main() {

	const double tau = 0.4;
	//const double sigma = 0.20; - WHY ARE THERE TWO SIGMAS?
	double meas[10000][2];
	
	read_data(meas);
	pdf(tau, meas);
	nll_tau(meas);

	return 0;
}	
