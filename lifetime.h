
//read in data
void read_data(double meas[][2]);
	
//fit function (p.d.f.) 
//- allows exploration of pdf for different taus and sigmas
void pdf(const double tau, const double sigma, const double meas[][2]); 

//function to get P for given tau and measurement
double get_P(const double tau, const double t, const double sigma);

//function to output nll for different tau values
void nll_tau(const double meas[][2]);

//function to find nll for given tau, t, sigma
double get_nll(const double tau, const double meas[][2]);

//initialise x and y values
void init(double x[], double y[], const double meas[][2]);

//initialise values for x and y arrays for cosh(x)
void init_cosh(double x[], double y[]);
	
//find coefficients for parabola eqn
void find_coeffs(double &A, double &B, const double x[], const double y[]);

//get location of parabolas minimum
void get_min(const double A, const double B, double x[], double y[], 
		const double meas[][2]);

//parabolic minimiser
void parabolic_minimiser(const double meas[][2]);

//find maximum y-value
double find_max(const double y[]);
	
//discard maximum y-value
void discard_max(double x[], double y[]);


//bisection method to find tau_plus and tau_minus from NLL
double bisect(const double nll_des, const double tau_outer, const double tau_inner, 
		const double meas[][2]);

//find standard deviations based on tau_plus & tau_minus
void get_stdev(double stdev[], const double xmin, const double y, 
		const double meas[][2]);

