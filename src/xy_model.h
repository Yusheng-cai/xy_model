#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace xy_model{
	double rotatespin(const double&);

	double randnum(double, double);

	std::vector<std::vector<double>> Random_initialize(int);

	double E_initial(const std::vector<std::vector<double>>&);

	double calc_deltaE(const std::vector<std::vector<double>>&,double,int,int);

	void run(int,int,double,std::string ENERGY="energy.dat",std::string CONFIG="config.dat",int printevery=1000);
}

