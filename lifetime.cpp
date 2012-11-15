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
		cout << "ERROR: read_data could not open datafile.txt \n";
		exit(1);
	}
	
	//2d array "meas" to contain time and error for each measurement
	for (int i = 0; i < 10000; i++) {
	      for(int j = 0; j < 2; j++) {
		      datafile >> meas[i][j];
	      }
	}
	cout << meas[1][1]; //TEST	
	datafile.close();
}

//calculate prob distribution P and output to file
void pdf(const double tau, double meas[][2]) {
	ofstream fitfile;
	fitfile.open("fitfunction.txt");

	//only proceed if fitfile opened successfully
	if (fitfile.fail()) {
		cout << "ERROR: fit could not open fitfunction.txt \n";
		exit(1);
	}

	//find P for each measurement for a given tau and output into "fitfunction.txt"
	for (int k = 0; k < 10000; k++) {
		double t = meas[k][0];
		double sigma = meas[k][1];
		
		fitfile << t << " " << sigma << " " << get_P << endl;
	}

	fitfile.close();

}

//function to find P for a given measurement
double get_P(double tau, double t, double sigma) {
	double err_input = ((sigma / tau) - (t / sigma)) / sqrt(2);
	double P = exp((sigma * sigma)/(2 * tau * tau) - (t / tau)) * erfc(err_input) / (2 * tau);

	return P;
}

//function to find NLL as a function of tau
void nll(double meas[][2]) {
	
