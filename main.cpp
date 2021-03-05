#include <iostream>
#include "simulation.h"
#include <ctime>
#include <stdlib.h>

int main(int argc, char** argv){
	std::srand(std::time(0));
	int N = atoi(argv[1]);
	double J = 0.2;

	run(N,0,J);


	
	

	return 0;
}
