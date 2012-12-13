#include "lifetime.h"
#include <fstream>
#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include <math.h> //to test for nans
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
	//specify range of t-values in order to integrate pdf
	//(want to be effectively integrating over all possible time)
	double max_t = 1000;
	double min_t = -10;
	double area = calculate_area(tau, sigma, max_t, min_t);

	//find P for each measurement for given tau and sigma 
	//normalise and output into "fit_function.txt"
	for (int i = 0; i < 10000; i++) {
		double t = measurements[i][0];
		double P = get_P_signal(tau, t, sigma); 

		double P_norm = P / area;
		fit_file << t << " " << P_norm << endl;
	}
	fit_file.close();
}

//find area
double calculate_area(const double tau, const double sigma, 
		const double upper_limit_t, const double lower_limit_t) {
	double A = evaluate_integral(tau, sigma, upper_limit_t) - 
		evaluate_integral(tau, sigma, lower_limit_t);
	cout << "Integral of PDF =  " << A << endl;
	return A;
}

//evaluate integral
double evaluate_integral(const double tau, const double sigma, 
		const double t) {
	double A = erf(t / (sqrt(2) * sigma)) / 2 -
		sqrt(exp((sigma * sigma / tau * tau) - 
					(2 * t / tau ))) * 
		erfc(((sigma / tau) - (t / sigma)) / sqrt(2));

	return A;
}

//find P for signal for a given measurement
double get_P_signal(const double tau, const double t, const double sigma) {
	const double err_input = ((sigma / tau) - (t / sigma)) / sqrt(2);
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

	//output warning if pdf for signal is negative
	if (abs(P_signal) != P_signal)
		cout << "WARNING: pdf for signal negative" << endl;

	P_background = get_P_background(t, sigma);
	P_total= a * P_signal + (1 - a) * P_background;

	//output warning if total pdf is negative
	if (abs(P_total) != P_total) { 
		cout << "WARNING: Total pdf negative" << endl;
		cout << "a = " << a << " tau = " << tau << endl;
	}
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
	for(int i = 0; i < i_max; i++) {//run through tau values	
		nll_file << tau_i << " " << get_nll(tau_i, measurements) << endl;	
		tau_i += d_tau;
	}	//finish running through tau values
	nll_file.close();
}

//find value of NLL for given tau
double get_nll(const double tau, const double measurements[][2]) {
	double nll = 0.00;
	//simulate having different sample sizes by changing "sample_size"
	const int sample_size = 10000;
	for(int i = 0; i < sample_size; i++) { //run through measurements
		int index = i % 10000;
		double t = measurements[index][0];
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
	//init_cosh(x, y);

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

	stdev_curvature(B);
	stdev_nll(xmin, y[3], measurements);
}

//find standard deivation based on curvature of parabolic estimate
void stdev_curvature(const double B) {
	cout << "\nSTANDARD DEVIATION BASED ON CURVATURE OF PARABOLA" << endl;
	cout << 1/sqrt(B) << endl;
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
	x[1] = 0.35;
	x[2] = 0.63;
	x[3] = 5.000;

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
	x[3] = (B * (x[0] + x[1]) - A) / (2. * B);
	y[3] = get_nll(x[3], measurements);
	//y[3] = cosh(x[3]);
}

//find value of maximum
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

	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	//put measured data into "par"
	double par[20000];
	measurements_to_par(par, measurements);

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	// Initialize method and iterate
	my_func.n = 2;
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = par;

	// Starting point
	x = gsl_vector_alloc(2);
	gsl_vector_set(x, 0, 0.2);
	gsl_vector_set(x, 1, 0.5);

	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc(T, 2);

	gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.0001, 0.1);

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 
				0.01);

		if (status == GSL_SUCCESS) {
			printf ("Minimum found at:\n");
			double a_mle = gsl_vector_get(s->x, 0);
			double tau_mle = gsl_vector_get(s->x, 1);

			//find standard deviations
			cout << "\nSTANDARD DEVIATIONS" << endl;
			double h = 0.00005;
			double dfda_da = get_dfda_da(h, a_mle, tau_mle, 
					par);
			double dfda_dtau = get_dfda_dtau(h, a_mle, 
					tau_mle, par);
			double dfdtau_dtau = get_dfdtau_dtau(h, a_mle, 
					tau_mle, par);
			double dfdtau_da = get_dfdtau_da(h, a_mle, 
					tau_mle, par);
			double det_H = dfda_da * 
			       dfdtau_dtau - dfda_dtau * dfda_dtau; 
			cout << "determinant = " << det_H << endl;
			double H11_inverse = dfdtau_dtau / det_H;
			double H12_inverse = - dfda_dtau / det_H;
			double H21_inverse = - dfda_dtau / det_H;
			double H22_inverse = dfda_da / det_H;

			double stdev_a = sqrt(H11_inverse);
			double stdev_tau = sqrt(H22_inverse);

			cout << "stdev_a = " << stdev_a << endl;
			cout << "stdev_tau = " << stdev_tau << endl;



		}
		/*printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		  iter, 
		  gsl_vector_get(s->x, 0),
		  gsl_vector_get(s->x, 1),
		  s->fval, size);*/
		cout << iter << " " << gsl_vector_get(s->x, 0) << 
			" " << 
			gsl_vector_get(s->x, 1) << " NLL=" << 
			s->f << endl;
		gsl_vector *grad = gsl_multimin_fdfminimizer_gradient(s);
		double g0 = gsl_vector_get(grad, 0);
		double g1 = gsl_vector_get(grad, 1);
		cout << "Gradient = " << g0 << ", " << g1 << 
			", " << endl;
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);
	return status;
}


//define function to be minimised
double my_f(const gsl_vector *v, void *params) {
	double a, tau;

	a = gsl_vector_get(v, 0);
	tau = gsl_vector_get(v, 1);
	double nll_total = get_nll_total(a, tau, params);	

	return nll_total;
}

//get nll total for case where it depends on tau and a
double get_nll_total(double a, const double tau, 
		void *params) {
	double *p = (double *)params;

	double nll_total = 0;
	int index_sigma = 0;
	int index_t = 0;
	for(int i = 0; i < 10000; i++) {
		index_t = 2 * i;
		index_sigma = 2 * i + 1;

		double t = p[index_t];
		double sigma = p[index_sigma];
		//a should not be greater than one
		//if it is, P_total becomes negative
		//and nll_total is nan
		if (a > 1)
			a -= 0.5;

		double P_total = get_P_total(a, tau, t, sigma);
		//TEST REAL NUMBER
		if (isnan(log(P_total))) {
			cout << "WARNING: LOG(PDF) IS NAN" << endl;
			cout << "a = " << a << " tau = " << tau << endl;
		}
		nll_total -= log(P_total);
	}
	return nll_total;
}

//gradient of f, df = (df/da, df/dtau)
void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
	double a, tau; 
	double h = 0.000001;

	a = gsl_vector_get(v, 0);
	tau = gsl_vector_get(v, 1);

	gsl_vector_set(df, 0, get_dfda(h, a, tau, params));
	gsl_vector_set(df, 1, get_dfdtau(h, a, tau, params));
}

//compute both fdf and df together
void my_fdf(const gsl_vector *x, void *params, double *f, 
		gsl_vector *df) {
	*f = my_f(x, params);
	my_df(x, params, df);
}


//put measured values into 1-D array "par"
void measurements_to_par(double par[], const double measurements[][2]) {
	for(int i = 0; i < 10000; i++) {
		int index_t = 2 * i;
		int index_sigma = 2 * i + 1;

		par[index_t] = measurements[i][0];
		par[index_sigma] = measurements[i][1];
	}
}

//use central difference scheme to estimate derivative of function
double get_dfda(const double h, const double a, const double tau, 
		void *params) {
	double dfda = (get_nll_total(a + h, tau, params) - 
			get_nll_total(a - h, tau, params)) / 
		(2 * h);
	return dfda;
}

double get_dfdtau(const double h, const double a, const double tau, 
		void *params) {
	double dfdtau = (get_nll_total(a, tau + h, params) - 
			get_nll_total(a, tau - h, params)) / 
		(2 * h);
	return dfdtau;
}

//approximate dfda_da using cds
double get_dfda_da(const double h, const double a, const double tau, 
		void *params) {
	double dfda_da = (get_dfda(h, a + h, tau, params) - 
			get_dfda(h, a - h, tau, params)) / (2 * h);
	return dfda_da;
}

//dfda_dtau
double get_dfda_dtau(const double h, const double a, const double tau, 
		void *params) {
	double dfda_dtau = (get_dfda(h, a, tau + h, params) - 
			get_dfda(h, a, tau - h, params)) / ( 2 * h);
	return dfda_dtau;
}

//dfdtau_dtau
double get_dfdtau_dtau(const double h, const double a, const double tau, 
		void *params) {
	double dfdtau_dtau = (get_dfdtau(h, a, tau + h, params) - 
			get_dfdtau(h, a, tau - h, params)) / ( 2 * h);
	return dfdtau_dtau;
}

//dfdtau_da
double get_dfdtau_da(const double h, const double a, const double tau, 
		void *params) {
	double dfdtau_dtau = (get_dfdtau(h, a + h, tau, params) - 
			get_dfdtau(h, a + h, tau, params)) / ( 2 * h);
	return dfdtau_dtau;
}
