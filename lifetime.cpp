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
//	cout << meas[2][1]; //TEST	
	datafile.close();
}

//calculate prob distribution P and output to file
void pdf(const double tau, const double sigma, const double meas[][2]) {
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
//		double sigma = meas[k][1]; - want to see effect of different sigmas
		
		fitfile << t << " " << get_P_sig(tau, t, sigma) << endl;
	}

	fitfile.close();

}

//function to find P (of signal) for a given measurement
double get_P_sig(const double tau, const double t, const double sigma) {
	double err_input = ((sigma / tau) - (t / sigma)) / sqrt(2);
	double P_sig = exp((sigma * sigma)/(2 * tau * tau) - 
			(t / tau)) * erfc(err_input) / (2 * tau);

	return P_sig;
}

//find P for background for a given measurement
double get_P_bkg(const double t, const double sigma) {
	const double pi = atan(1) *4;
	double P_bkg = (exp((-t * t) / (2 * sigma * sigma))) / (sigma * sqrt(2 * pi));

	return P_bkg;
}

//find P for background and signal for a given measurement
double get_P_total(const double a, const double tau, const double t, const double sigma) {
	double P_sig, P_bkg, P_tot;
	P_sig = get_P_sig(tau, t, sigma);
	P_bkg = get_P_bkg(t, sigma);
	P_tot = a * P_sig + (1 - a) * P_bkg;
	
	return P_tot;
}


//function to output NLL for different tau values
void nll_tau(const double meas[][2]) {
	double tau_min = 0.42;
	double d_tau = 0.00001;
	double tau_max = 0.431;
	//find appropriate k_max for desired tau_max
	double k_max = (tau_max - tau_min) / d_tau; 
	//always round up so desired range is included 
	k_max = round(0.60 + k_max); 

	ofstream nllfile;
	nllfile.open("nllfunction.txt");

	//only proceed if fitfile opened successfully
	if (nllfile.fail()) {
		cout << "ERROR: find_nll could not open nllfunction.txt \n";
		exit(1);
	}

	double tau_k = tau_min;
	for(int k = 0; k < k_max; k++) { //run through tau values	
		nllfile << tau_k << " " << get_nll(tau_k, meas) << endl;	
		tau_k += d_tau;
	} //finish running through tau values
	nllfile.close();
}

//find value of NLL for given tau
double get_nll(const double tau, const double meas[][2]) {
	double nll = 0.00;
	int sample_size = 10000;
		for(int i = 0; i < sample_size; i++) { //run through measurements
			int index = i % 10000;
			double t = abs(meas[index][0]);
			double sigma = meas[index][1];

			double P_sig = get_P_sig(tau, t, sigma);
			nll -= log(P_sig);
		} //finish running through measurements
	return nll;
}	

//parabolic minimiser
void parabolic_minimiser(const double meas[][2]) {
	cout << "\nPARABOLIC MINIMISATION" << endl;
	double x[4];
	double y[4];
	
	double A = 0;
	double B = 0;
	
	//count number of iterations
	int iterations =0;
	double xmin;
	double xmin_prev;

	//initialise x and y values
	init(x, y, meas);

	do {
		find_coeffs(A, B, x, y);
		xmin_prev = x[3];
		get_min(A, B, x, y, meas);
		xmin = x[3]; 
		discard_max(x, y);
		iterations++;
	} 
	//specify convergence criterion
	while (abs(xmin - xmin_prev) > 0.000001);

	cout << "Minimum value of NLL = " << y[3] << endl;
	cout << "tau-value at minimum = " << xmin << endl;
	cout << "Number of iterations = " << iterations << endl;

	stdev_parabolic(A, B, xmin, x[0], x[1], y[0], y[3]);
	stdev_nll(xmin, y[3], meas);
}

//find and output standard deviation based on latest parabolic estimate
void stdev_parabolic(const double A, const double B, const double xmin, 
		const double x0, const double x1, const double y0, 
		const double y3) {
	cout << "\nSTANDARD DEVIATION BASED ON LATEST PARABOLIC ESTIMATE" << endl;
	double nll_stdv = y3 + 0.5;
	const double a = B;
	const double b = A - B * x0 - B * x1;
	const double c = B * x0 * x1 - A * x0 + y0 - nll_stdv;

	//use quadratic formula to find tau value at nll_stdv
	double stdev = (-b + sqrt(b * b - 4 * a * c)) / (2 * a) - xmin;

	cout << "The NLL-value being used is " << nll_stdv << endl;
	cout << "stdev_para = " << stdev << endl;
}

//find and output standard deviations based on tau_plus & tau_minus
void stdev_nll(const double xmin, const double y3, const double meas[][2]) {
	cout << "\nSTANDARD DEVIATION BASED ON CHANGE IN NLL" << endl;
	const double nll_stdev = y3 + 0.5;
	double stdev[2];
	
	//find tau_plus
	double tau_left = xmin;
	double tau_right = 5.5;
	double tau_plus = bisect(nll_stdev, tau_left, tau_right, meas);
	cout << "tau_plus = " << tau_plus << "\n" << endl;
	
	//find tau_minus
	tau_left = 0.05;
	tau_right = xmin;
	double tau_minus = bisect(nll_stdev, tau_left, tau_right, meas);
	cout << "tau_minus = " << tau_minus << "\n" << endl;

	double stdev_plus = tau_plus - xmin;
	double stdev_minus = xmin - tau_minus;
	stdev[0] = stdev_plus;
	stdev[1] = stdev_minus;

	cout << "stdev_plus = " << stdev[0] << endl;
	cout << "stdev_minus = " << stdev[1] << endl;

}

//bisection method to find tau_plus and tau_minus from NLL
double bisect(const double nll_des, double tau_left, double tau_right, 
		const double meas[][2]) {
	double tau_mid;
	double nll_mid;
	const double tolerance = 0.001;
	double i_max = 40;
	int evals;
	
	double nll_left = get_nll(tau_left, meas);
	double nll_right = get_nll(tau_right, meas);
	
	for (int i = 0; i < i_max; i++) {
		tau_mid = (tau_left + tau_right) / 2;
		nll_mid = get_nll(tau_mid, meas);

		//check if nll is within tolerance of desired nll value
		if (abs(nll_mid - nll_des) <= tolerance)
			break;
		//if finding tau_minus
		else if (nll_left > nll_right) {
			if (nll_mid < nll_des)
				tau_right = tau_mid;
			else tau_left = tau_mid;
		}
		//if finding tau_plus
		else {
			if (nll_mid > nll_des)
				tau_right = tau_mid;
			else tau_left = tau_mid;
		}

		// quit & error message if large number of iterations but no solution
		if (i == i_max - 1) {
			cout << "No solution after " << i_max << " iterations" << endl;
			cout << "NLL_mid = " << nll_mid << endl;
			cout << "tau_mid = " << tau_mid << endl;
			return 100;
		}
	evals = i;
	}
	//check nll value
	cout << "Function evaluated " << evals + 1 << " times." << endl;
	cout << "The NLL-value settled on is " << nll_mid << endl;
	return tau_mid;
}

//initialise values for x and y arrays
void init(double x[], double y[], const double meas[][2]) {
	x[0] = 0.22;
	x[1] = 0.55;
	x[2] = 0.63;
	x[3] = 6.000;
	
	y[0] = get_nll(x[0], meas);
	y[1] = get_nll(x[1], meas);
	y[2] = get_nll(x[2], meas);
	y[3] = get_nll(x[3], meas);
}

//initialise values for x and y arrays for cosh(x)
void init_cosh(double x[], double y[]) {
	x[0] = -1.5;
	x[1] = 0.5;
	x[2] = 1.0;
	x[3] = 100;

	y[0] = cosh(x[0]);
	y[1] = cosh(x[1]);
	y[2] = cosh(x[2]);
	y[3] = cosh(x[3]);
}

//find coefficients for the quadratic function which connects x0, x1, x2
void find_coeffs(double &A, double &B, const double x[], const double y[]) {
	A = (y[1] - y[0]) / (x[1] - x[0]);
	B = ((y[2] - y[0]) / (x[2] - x[0]) - (y[1] - y[0]) / (x[1] - x[0])) / (x[2] - x[1]);
}
	

//get location of the parabola's minimum
void get_min(const double A, const double B, double x[], double y[], const double meas[][2]) {
	x[3] = (B * (x[0] + x[1]) - A) / (2 * B);
	y[3] = get_nll(x[3], meas);
	//y[3] = cosh(x[3]);
}

//find location of maximum
double find_max(const double y[]) {
	double y_max = y[0];
	if (y[1] > y_max)
		y_max = y[1];
	if (y[2] > y_max)
		y_max = y[2];
	if (y[3] > y_max)
		y_max = y[3];
	return y_max;
}

//discard maximum y-value
void discard_max(double x[], double y[]) {
	double y_old[4];
	for (int i = 0; i < 4; i++)
		y_old[i] = y[i];

	double x_old[4];
	for (int i = 0; i < 4; i++)
		x_old[i] = x[i];


	double y_max = find_max(y_old);
	if (y_max == y_old[0]) {
		y[0] = y_old[1];
		y[1] = y_old[2];
		y[2] = y_old[3];

		x[0] = x_old[1];
		x[1] = x_old[2];
		x[2] = x_old[3];
	}
	else if (y_max == y_old[1]) {
		y[0] = y_old[0];
		y[1] = y_old[2];
		y[2] = y_old[3];

		x[0] = x_old[0];
		x[1] = x_old[2];
		x[2] = x_old[3];
	}
	else if (y_max == y_old[2]) {
		y[0] = y_old[0];
		y[1] = y_old[1];
		y[2] = y_old[3];

		x[0] = x_old[0];
		x[1] = x_old[1];
		x[2] = x_old[3];
	}
	else if (y_max == y_old[3]) {
		y[0] = y_old[0];
		y[1] = y_old[1];
		y[2] = y_old[2];

		x[0] = x_old[0];
		x[1] = x_old[1];
		x[2] = x_old[2];
	}
}

//MULTIDIMENSIONAL MINIMISER
void multimin() {
	gsl_vector *v = gsl_vector_alloc(20000);
}
//read measurements into a gsl_vector v
//elements alternate between t and sigma
void meas_to_vector() {
	gsl_vector *v = gsl_vector_alloc(20000);
	{
	FILE *f = fopen("lifetime.txt", "r");
	gsl_vector_fscanf(f, v);	
	fclose(f);
	}
	for(int i = 0; i < 10; i++) {
		printf("%g\n", gsl_vector_get(v,i));
	}
	gsl_vector_free(v);
}
