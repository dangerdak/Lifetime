#include "lifetime.h"
#include <fstream>
#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include <cstdlib> //so I can use "exit"
#include <iostream>

using namespace std;

//read in data to 2-d array
void read_data(double meas[][2]) {
	ifstream datafile;
	datafile.open("lifetime.txt");

	//only proceed if datafile opened successfully
	if (datafile.fail()) {
		cout << "ERROR: read_data could not open lifetime.txt \n";
		exit(1);
	}
	
	//2d array "meas" to contain time and error for each measurement
	for (int i = 0; i < 10000; i++) {
	      for(int j = 0; j < 2; j++) {
		      datafile >> meas[i][j];
	      }
	}
	cout << meas[2][1]; //TEST	
	datafile.close();
}

//calculate prob distribution P and output to file
void pdf(const double tau, const double meas[][2]) {
	ofstream fitfile;
	fitfile.open("fitfunction.txt");

	//only proceed if fitfile opened successfully
	if (fitfile.fail()) {
		cout << "ERROR: fit could not open fitfunction.txt \n";
		exit(1);
	}

	//find P for each measurement for a given tau and output into "fitfunction.txt"
	for (int k = 0; k < 10000; k++) {
		double t = abs(meas[k][0]);
		double sigma = meas[k][1];
		
		fitfile << t << " " << sigma << " " << get_P(tau, t, sigma) << endl;
	}

	fitfile.close();

}

//function to find P for a given measurement
double get_P(const double tau, const double t, const double sigma) {
	double err_input = ((sigma / tau) - (t / sigma)) / sqrt(2);
	double P = exp((sigma * sigma)/(2 * tau * tau) - (t / tau)) * erfc(err_input) / (2 * tau);

	return P;
}

//function to output NLL for different tau values
void nll_tau(const double meas[][2]) {
	double tau_min = 0.42;
	double d_tau = 0.00001;
	double tau_max = 0.44;
	double k_max = (tau_max - tau_min) / d_tau; //find appropriate k_max for desired tau_max
	k_max = round(0.60 + k_max); //always round up so desired range is included 

	ofstream nllfile;
	nllfile.open("nllfunction.txt");

	//only proceed if fitfile opened successfully
	if (nllfile.fail()) {
		cout << "ERROR: find_nll could not open nllfunction.txt \n";
		exit(1);
	}

	double tau_k = tau_min;
	for(int k = 0; k < k_max; k++) { //run through tau values	
		double nll = 0.00;
		for(int i = 0; i < 10000; i++) { //run through measurements
			double t = abs(meas[i][0]);
			double sigma = meas[i][1];

			double P = get_P(tau_k, t, sigma);
			nll -= log(P);
		} //finish running through measurements
		nllfile << tau_k << " " << nll << endl;	
		tau_k += d_tau;
	} //finish running through tau values
	nllfile.close();
}

//find value of NLL for given tau, t, sigma
double get_nll(const double tau, const double t, const double sigma) {
