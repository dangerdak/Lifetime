#include "lifetime.h"

using namespace std;

int main() {

	//constants used with calculate_pdf
	const double tau = 0.40;
	const double sigma = 0.20;
	
	//constants used with nll_vs_tau
	//tau_min and tau_max must bracket minimum
	const double tau_min = 0.400;
	const double d_tau = 0.0001;
	const double tau_max = 0.410;

	double measurements[10000][2];
	
	read_data(measurements);
	calculate_pdf(tau, sigma, measurements);
	nll_vs_tau(tau_min, d_tau, tau_max, measurements);
//	parabolic_minimiser(measurements);
//	multimin(measurements);
	return 0;
}	
