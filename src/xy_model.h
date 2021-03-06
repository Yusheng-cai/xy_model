#include <vector>

namespace xy_model{
	double rotatespin(const double &spin);
	double randnum(float low, float high);
	std::vector<std::vector<double>> Random_initialize(int N);
	double E_initial(const std::vector<std::vector<double>> &lattice,double J);
	double calc_deltaE(const std::vector<std::vector<double>> &lattice,double updated_spin,int i,int j,double J);
	void run(int N,int nsweeps,double J,double T);
}

