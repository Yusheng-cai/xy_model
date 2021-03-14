#include <iostream>
#include <ctime>
#include <stdlib.h>
#include "xy_model.h"
#include <string>

int main(int argc, char** argv){
	std::srand(std::time(0));
	int N = atoi(argv[1]);
	int nsweeps = atoi(argv[2]);
	double T = atof(argv[3]);
	double J = atof(argv[4]);

	if(argc > 5){	
		std::string energyname = argv[5];
		std::string configname = argv[6];
		xy_model::run(N,nsweeps,J,T,energyname,configname);
	}
	else{
		xy_model::run(N,nsweeps,J,T);
	}


	return 0;
}
