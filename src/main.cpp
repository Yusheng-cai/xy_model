#include <iostream>
#include <ctime>
#include <stdlib.h>
#include "xy_model.h"

int main(int argc, char** argv){
	std::srand(std::time(0));
	int N = atoi(argv[1]);
	double J = 2;
	double T = 300;

	xy_model::run(N,0,J,T);
	return 0;
}
