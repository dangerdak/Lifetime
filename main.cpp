#include "lifetime.h"

using namespace std;

int main() {

	//const double tau = 0.01;
	//const double sigma = 0.2;
	double measurements[10000][2];
	
	read_data(measurements);
	//calculate_pdf(tau, sigma, measurements);
	//nll_tau(measurements);
	parabolic_minimiser(measurements);
	multimin(measurements);
	return 0;
}	
