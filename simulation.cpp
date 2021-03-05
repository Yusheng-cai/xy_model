#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#define PI 3.14159265


double randomangle(){
	// Generate a random angle between -180 and 180
	double RANDANGLE = -180.0 + (double) 360.0*rand()/RAND_MAX;
	
	return RANDANGLE;
}

double rotatespin(const double &angle){
	double rotated_angle;
	double RANDANGLE = randomangle();
	
	//Rotate clockwise if it the negative 
	if ((angle + RANDANGLE) < -180){
		rotated_angle = angle + RANDANGLE + 360.0;
	}
	else if ((angle + RANDANGLE) > 180){
		rotated_angle = angle + RANDANGLE - 360.0;	
	}
	else{
		rotated_angle = angle +RANDANGLE;
	}

	return rotated_angle;
}


double E_initial(std::vector<std::vector<double>> &lattice,double J){
	/*
	 * Function that calculates the energy of the initial spin system 
	 *
	 * Args:
	 * lattice(std::vector<std::vector<double>>): lattice of the angles on the lattice formed with positive x direction
	 * J(double): The energy constant J*cos(thetai-thetaj).
	 * 
	 *
	 * Returns:
	 * Energy(double): The energy as a double
	 */
	int rows,cols;
	double lval,rval,dval,uval,currval;
	double energy=0.0;

	// Get the sizes of the vector<vector<double>>
	rows = lattice.size();
	cols = lattice[0].size();

	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			currval = lattice[i][j];
			// Check the upper pbc
			if(i == 0){
				uval=lattice[rows-1][j];
			}	
			else{
				uval=lattice[i-1][j];
			}

			// Check the lower pbc
			if(i == rows-1){
				dval=lattice[0][j];
			}
			else{
				dval=lattice[i+1][j];
			}

			// Check the left pbc
			if(j == 0){
				lval=lattice[i][cols-1];
			}
			else{
				lval=lattice[i][j-1];
			}

			// Check the right pbc
			if(j == cols-1){
				rval = lattice[i][0];
			}
			else{
				rval = lattice[i][j+1];
			}

			energy += -J*(cos((currval-rval)*PI/180.0)+cos((currval-lval)*PI/180.0)+cos((currval-uval)*PI/180.0)+\
					cos((currval-dval)*PI/180.0));
		}
	}

	return energy;
}

void run(int N,int nsweeps,double J){
	std::vector<std::vector<double>> lattice(N,std::vector<double>(N,0));

	// Randomly initialize lattice
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			lattice[i][j] = randomangle();
			std::cout << lattice[i][j] << std::endl;
		}
	}

	double E = E_initial(lattice,J);
	std::cout << E << std::endl;
}
