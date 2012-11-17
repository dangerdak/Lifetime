
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

//find coefficients for parabola eqn
void find_coeffs(double &A, double &B, const double x0, const double x1, 
	const double x2, const double y0, const double y1, const double y2); 
	
