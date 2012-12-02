#include "lifetime.h"
#include <fstream>
#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include <cstdlib> //so I can use "exit"
#include <iostream>
#include <gsl/gsl_multimin.h>

using namespace std;

//read in data to 2-d array
void read_data(double measurements[][2]) {
	ifstream data_file;
	data_file.open("lifetime.txt");

	//only proceed if data file opened successfully
	if (data_file.fail()) {
		cerr << "ERROR: read_data could not open lifetime.txt \n";
		exit(1);
	}
	//put time measurements in column 0
	//and associated errors in column 1
	for (int i = 0; i < 10000; i++) 
	      for(int j = 0; j < 2; j++) 
		      data_file >> measurements[i][j];
	      
	data_file.close();
}

//calculate probability distribution P for given tau and sigma, and output to file
void calculate_pdf(const double tau, const double sigma, 
		const double measurements[][2]) {
	ofstream fit_file;
	fit_file.open("fit_function.txt");

	//only proceed if fit file opened successfully
	if (fit_file.fail()) {
		cerr << "ERROR: calculate_pdf could not open fit_function.txt \n";
		exit(1);
	}
	//find P for each measurement for given tau and sigma 
	//and output into "fit_function.txt"
	for (int k = 0; k < 10000; k++) {
		double t = abs(measurements[k][0]);
		fit_file << t << " " << get_P_signal(tau, t, sigma) << endl;
	}
	fit_file.close();
}

//find P for signal for a given measurement
double get_P_signal(const double tau, const double t, const double sigma) {
	double err_input = ((sigma / tau) - (t / sigma)) / sqrt(2);
	double P_signal = exp((sigma * sigma)/(2 * tau * tau) - 
			(t / tau)) * erfc(err_input) / (2 * tau);

	return P_signal;
}

//find P for background for a given measurement
//M_PI is value of pi defined in gsl
double get_P_background(const double t, const double sigma) {
	double P_background = (exp((-t * t) / (2 * sigma * sigma))) / 
		(sigma * sqrt(2 * M_PI));

	return P_background;
}

//find total P for background and signal for a given measurement
double get_P_total(const double a, const double tau, const double t, const double sigma) {
	double P_signal, P_background, P_total;
	P_signal = get_P_signal(tau, t, sigma);
	P_background = get_P_background(t, sigma);
	P_total= a * P_signal + (1 - a) * P_background;
	
	return P_total;
}


//function to output NLL for different tau values
void nll_vs_tau(const double tau_min, const double d_tau, 
		const double tau_max, const double measurements[][2]) {
	//find appropriate number of iterations for desired tau_max
	//always round up so desired range is included 
	const int i_max = ceil((tau_max - tau_min) / d_tau); 

	ofstream nll_file;
	nll_file.open("nll_function.txt");

	//only proceed if nll_file opened successfully
	if (nll_file.fail()) {
		cerr << "ERROR: find_nll could not open nll_function.txt \n";
		exit(1);
	}

	double tau_i = tau_min;
	for(int i = 0; i < i_max; i++) //run through tau values	
		nll_file << tau_i << " " << get_nll(tau_i, measurements) << endl;	
	//finish running through tau values
	nll_file.close();
}

//find value of NLL for given tau
double get_nll(const double tau, const double measurements[][2]) {
	double nll = 0.00;
	//simulate having different sample sizes by changing "sample_size"
	const int sample_size = 10000;
		for(int i = 0; i < sample_size; i++) { //run through measurements
			int index = i % 10000;
			double t = abs(measurements[index][0]);
			double sigma = measurements[index][1];

			double P_signal = get_P_signal(tau, t, sigma);
			nll -= log(P_signal);
		} //finish running through measurements
	return nll;
}	

//parabolic minimiser
//using parabola described by
//g(x) = y[0] + A*(x - x[0]) + B*(x - x[0])*(x - x[1])
void parabolic_minimiser(const double measurements[][2]) {
	cout << "\nPARABOLIC MINIMISATION" << endl;
	//arrays x and y contain coordinates 
	//for the 3 points to which the parabola is fitted
	//and for the minimum of the parabola (in the final element)
	double x[4];
	double y[4];
	
	double A = 0;
	double B = 0;
	
	//count number of iterations
	int iterations = 0;
	double xmin;
	double xmin_previous;

	//initialise x and y values
	init(x, y, measurements);

	do {
		find_coeffs(A, B, x, y);
		xmin_previous = x[3];
		
		get_min(A, B, x, y, measurements);
		xmin = x[3]; 

		discard_max(x, y);
		iterations++;
	} 
	//specify convergence criterion
	while (abs(xmin - xmin_previous) > 0.000001);

	cout << "Minimum value of NLL = " << y[3] << endl;
	cout << "tau-value at minimum = " << xmin << endl;
	cout << "Number of iterations = " << iterations << endl;

	stdev_parabolic(A, B, xmin, x[0], x[1], y[0], y[3]);
	stdev_nll(xmin, y[3], measurements);
}

//find and output standard deviation based on latest parabolic estimate
void stdev_parabolic(const double A, const double B, const double xmin, 
		const double x0, const double x1, const double y0, 
		const double y3) {
	cout << "\nSTANDARD DEVIATION BASED ON LATEST PARABOLIC ESTIMATE" << endl;
	const double nll_stdv = y3 + 0.5;

	//use quadratic formula to find tau value at nll_stdv
	//for function of form ax^2 + bx + c = 0
	const double a = B;
	const double b = A - B * x0 - B * x1;
	const double c = B * x0 * x1 - A * x0 + y0 - nll_stdv;
	const double stdev = (-b + sqrt(b * b - 4 * a * c)) / (2 * a) - xmin;

	cout << "The NLL-value being used is " << nll_stdv << endl;
	cout << "stdev_parabolic_estimate = " << stdev << endl;
}

//find and output standard deviations based on change in NLL
void stdev_nll(const double xmin, const double y3, const double measurements[][2]) {
	cout << "\nSTANDARD DEVIATION BASED ON CHANGE IN NLL" << endl;
	const double nll_stdev = y3 + 0.5;
	double stdev[2];
	
	//find tau_plus
	double tau_left = xmin;
	double tau_right = 5.5;
	double tau_plus = bisect(nll_stdev, tau_left, tau_right, measurements);
	cout << "tau_plus = " << tau_plus << "\n" << endl;
	
	//find tau_minus
	tau_left = 0.05;
	tau_right = xmin;
	double tau_minus = bisect(nll_stdev, tau_left, tau_right, measurements);
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
		const double measurements[][2]) {
	double tau_mid;
	double nll_mid;
	const double tolerance = 0.001;
	double i_max = 40;
	int evals;
	
	double nll_left = get_nll(tau_left, measurements);
	double nll_right = get_nll(tau_right, measurements);
	
	for (int i = 0; i < i_max; i++) {
		tau_mid = (tau_left + tau_right) / 2;
		nll_mid = get_nll(tau_mid, measurements);

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
void init(double x[], double y[], const double measurements[][2]) {
	x[0] = 0.22;
	x[1] = 0.55;
	x[2] = 0.63;
	x[3] = 6.000;
	
	y[0] = get_nll(x[0], measurements);
	y[1] = get_nll(x[1], measurements);
	y[2] = get_nll(x[2], measurements);
	y[3] = get_nll(x[3], measurements);
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
void get_min(const double A, const double B, double x[], double y[], const double measurements[][2]) {
	x[3] = (B * (x[0] + x[1]) - A) / (2 * B);
	y[3] = get_nll(x[3], measurements);
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
int multimin(const double measurements[][2]) {
	cout << "\n2-D MINIMISATION" << endl;
	//put measured data into "par"
	double par[20000];
	measurements_to_par(par, measurements);

	const gsl_multimin_fminimizer_type *T = 
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;

	// Starting point - WHAT SHOULD THE X-VALUES BE????
	x = gsl_vector_alloc(2);
	gsl_vector_set(x, 0, 0.7);
	gsl_vector_set(x, 1, 0.5);

	// Set initial step sizes to 0.2
	ss = gsl_vector_alloc(2);
	gsl_vector_set_all(ss, 0.2);

	// Initialize method and iterate
	minex_func.n = 2;
	minex_func.f = my_f;
	minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc(T, 2);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
	
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 0.01);

		if (status == GSL_SUCCESS) {
			printf ("Converged to minimum at:\n");
		}
		/*printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
				iter, 
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				s->fval, size);*/
		cout << iter << " " << gsl_vector_get(s->x, 0) << " " << 
			gsl_vector_get(s->x, 1) << " f()=" << s->fval << " size=" << 
			size << endl;
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	return status;
}


//define function to be minimised
double my_f(const gsl_vector *v, void *params) {
	double a, tau;
	double *p = (double *)params;

	a = gsl_vector_get(v, 0);
	tau = gsl_vector_get(v, 1);

	double nll_total = 0;
	int index_sigma = 0;
	int index_t = 0;
	for(int i = 0; i < 10000; i++) {
		index_t = 2 * i;
		index_sigma = 2 * i + 1;

		double t = abs(p[index_t]);
		double sigma = p[index_sigma];

		double P_total = get_P_total(a, tau, t, sigma);
		nll_total -= log(P_total);
	}

	return nll_total;
}

//put measured values into array "par"
void measurements_to_par(double par[], const double measurements[][2]) {
	for(int i = 0; i < 10000; i++) {
		int index_t = 2 * i;
		int index_sigma = 2 * i + 1;

		par[index_t] = measurements[i][0];
		par[index_sigma] = measurements[i][1];
	}
}
