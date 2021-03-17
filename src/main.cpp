#include "xy_model.h"
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <unistd.h>

int main(int argc, char* argv[]){
	std::srand(std::time(0));

	// command line arguments must be larger than 3 for N, nsweeps and kbTJ
	if(argc < 4){
		std::cerr << "You must provide: N(number of grids), nsweeps(number of sweeps to be performed), kbTJ(kbT/J)" << std::endl;
	       	return 1;	
	}
	
	int opt;
	int N=0;
	int nsweeps=0;
	int printevery=1000;
	double kbTJ=0.0;
	std::string energyname="energy.dat";
	std::string propertyname="property.dat";
	std::string configname="config.dat";

	
	for(;;){	
		switch(getopt(argc,argv,"n:s:k:")){
		case 'n':
			N=atoi(optarg);
			continue;
		case 's':
			nsweeps=atoi(optarg);
			continue;
		case 'k':
			kbTJ=atof(optarg);
			break;
			}
		break;
	}
	
	std::cout << N << std::endl;
	std::cout << nsweeps << std::endl;
	std::cout << kbTJ << std::endl;
	xy_model::run(N,nsweeps,kbTJ,energyname,configname,propertyname,printevery);
	

	return 0;
}
