#include <gsl/gsl_vector.h>

//read in data from "lifetime.txt" to array "measurements"
//time in column 0, associated errors in column 1
void read_data(double measurements[][2]);
	
//calculate probability distribution function P
//for given average time (tau) and error (sigma) for each measurement
//output data points to file "fit_function.txt"
//with time in column 0, value of pdf in column 1
void calculate_pdf(const double tau, const double sigma, 
		const double measurements[][2]); 
//find max and min t
void find_max_min_t(double &max_t, double &min_t, 
		const double measurements[][2]);

//find area
double calculate_area(const double tau, const double sigmai, 
		const double upper_limit_t, const double lower_limit_t);
	
//evaluate integral
double evaluate_integral(const double tau, const double sigma, 
		const double t);

//find P (of signal) for given tau and measurement
double get_P_signal(const double tau, const double t, const double sigma);

//find P for background for a given measurement
double get_P_background(const double t, const double sigma);

//find total P for background and signal combined for a given measurement
double get_P_total(const double a, const double tau, const double t, 
		const double sigma);
	
//output nll for different tau values
//tau_min and tau_max must bracket minimum
//d_tau is interval between tau-values
void nll_vs_tau(const double tau_min, const double d_tau, 
		const double tau_max, const double measurements[][2]);

//find value of NLL for given tau
double get_nll(const double tau, const double measurements[][2]);

//parabolic minimiser
//using parabola described by
//g(x) = y[0] + A*(x - x[0]) + B*(x - x[0])*(x - x[1])
void parabolic_minimiser(const double measurements[][2]);

//initialise x and y values for parabolic minimiser
void init(double x[], double y[], const double measurements[][2]);

//initialise values for x and y arrays to minimise cosh(x)
//using parabolic minimiser
void init_cosh(double x[], double y[]);
	
//find coefficients (A and B) for parabola equation
void find_coeffs(double &A, double &B, const double x[], 
		const double y[]);

//get location of parabolas minimum
void get_min(const double A, const double B, double x[], double y[], 
		const double measurements[][2]);

//find standard deivation based on curvature of parabolic estimate
void stdev_curvature(const double B);

//find maximum value in array y
double find_max(const double y[]);
	
//discard maximum y-value
void discard_max(double x[], double y[]);

//find and output standard deviations based on change in NLL
void stdev_nll(const double xmin, const double y3, const double measurements[][2]);

//bisection method to find tau_plus and tau_minus from NLL
//i.e. the tau-values one standard deviation above/below the minimum
double bisect(const double nll_des, const double tau_outer, const double tau_inner, 
		const double measurements[][2]);

//MULTIDIMENSIONAL MINIMISER
//minimise function described by "my_f"
int multimin(const double measurements[][2]);

//read measurements into a gsl_vector which alternates between t and sigma
void measurements_to_vector(gsl_vector *v);
	
//define function to be minimised
double my_f(const gsl_vector *v, void *params);

//get nll total for case where it depends on tau and a
double get_nll_total(double a, const double tau, 
		void *params);

//gradient of f, df = (df/da, df/dtau)
void my_df(const gsl_vector *v, void *params, gsl_vector *df);

//compute both fdf and df together
void my_fdf(const gsl_vector *x, void *params, double *f, 
		gsl_vector *df);

//put measured values into 1-D array "par"
void measurements_to_par(double par[], const double measurements[][2]);

//use central difference scheme to estimate derivative of function
double get_dfda(const double h, const double a, const double tau, 
		void *params);
	
double get_dfdtau(const double h, const double a, const double tau, 
		void *params);
	
