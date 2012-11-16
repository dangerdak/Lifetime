
//read in data
void read_data(double meas[][2]);
	
//fit function (p.d.f.)
void pdf(const double tau, const double meas[][2]); 

//function to get P for given tau and measurement
double get_P(const double tau, const double t, const double sigma);

//function to output nll for different tau values
void nll_tau(const double meas[][2]);
