#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <type_traits>
#include <vector>
#include <omp.h>

namespace xy_model{
	double rotatespin(const double&);

	double randnum(double, double);

	std::vector<std::vector<double>> Random_initialize(int);

	double E_initial(const std::vector<std::vector<double>>&);

	double calc_deltaE(const std::vector<std::vector<double>>&,double,int,int);

	void run(int,int,double,double,std::string,std::string,std::string,int printevery=1000);
}

