#include "xy_model.h"
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <string>

int main(int argc, char** argv){
	std::srand(std::time(0));
	int N = atoi(argv[1]);
	int nsweeps = atoi(argv[2]);
	double kbTJ = atof(argv[3]);

	if(argc > 4){	
		std::string energyname = argv[4];
		std::string configname = argv[5];
		xy_model::run(N,nsweeps,kbTJ,energyname,configname);
	}
	else{
		xy_model::run(N,nsweeps,kbTJ);
	}


	return 0;
}
