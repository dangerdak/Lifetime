
//read in data
void read_data(double meas[][2]);
	
//fit function (p.d.f.)
void pdf(const double tau, double meas[][2]); 

//function to get P for given tau and measurement
double get_P(double tau, double t, double sigma);

//function to output nll for different tau values
void find_nll(double meas[][2]);
