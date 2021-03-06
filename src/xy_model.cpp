#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "xy_model.h"
#include <omp.h>


#define PI 3.14159265
#define RAD 3.14159265/180.0


double xy_model::randomangle(){
	/*
	 * Function that generates a random angle between -180 and 180
	 *
	 * Returns:
	 * RANDANGLE(double): A random number between -180 and 180
	 */
	double RANDANGLE = -180.0 + (double) 360.0*rand()/RAND_MAX;
	
	return RANDANGLE;
}

double xy_model::rotatespin(const double &angle){
	/*
	 * Function that rotates an spin on the lattice by a random angle 
	 * will return a number between -180 and 180 for sure
	 *
	 * Args:
	 * angle(double): The original angle on the lattice
	 * 
	 * Return:
	 * rotated_angle(double): The rotated angle by a random number of degrees
	 */
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


std::vector<std::vector<double>> xy_model::Random_initialize(int N){
	/*
	 * Function that randomly initializes the configuration of the lattice at time 0
	 *
	 * Args:
	 * N(double): The number of discretized element
	 *
	 * Return:
	 * lattice(std::vector<std::vector<double>>): The lattice where i,j represents i,j element's angle with the positive x axis
	 */
	std::vector<std::vector<double>> lattice(N,std::vector<double>(N,0));

	// Randomly initialize lattice
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			lattice[i][j] = randomangle();
		}
	}

	return lattice;	
}

double xy_model::E_initial(const std::vector<std::vector<double>> &lattice,double J){
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

	double start=omp_get_wtime();

	#pragma omp parallel for shared(lattice) reduction(+:energy)	
	for(int k=0;k<rows*cols;k++){
		int i=k/rows;
		int j=k%rows;
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

		energy +=-J*(cos((currval-rval)*RAD)+cos((currval-lval)*RAD)+cos((currval-uval)*RAD)+\
				cos((currval-dval)*RAD));	
	}
	double end=omp_get_wtime();
	std::cout << end-start<<std::endl;
	


	return energy;
}

double xy_model::calc_deltaE(const std::vector<std::vector<double>> &lattice, double updated_spin,int i, int j,double J){
	/*
	 * Function that calculates the change of energy from changing the spin
	 * 
	 * Args:
	 * lattice(std::vector<std::vector<double>>): A lattice that holds all the spins
	 * updated_spin(double): The value of the spin to be updated
	 * i(int): The row that the spin belongs to
	 * j(int): The column that the spin belongs to
	 *
	 * Return:
	 * deltaE(double): The change of energy upon changing the spin
	 */
	double rows=lattice.size();
	double cols=lattice[0].size();
	double uval,dval,lval,rval;
	
	double currval = lattice[i][j];

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

	double dE =J*(cos((currval-rval)*RAD) - cos((updated_spin-rval)*RAD)+\
			cos((currval-lval)*RAD) - cos((updated_spin-lval)*RAD)+\
			cos((currval-uval)*RAD) - cos((updated_spin-uval)*RAD)+cos((currval-dval)*RAD) - cos((updated_spin - dval)*RAD));	
	
       return dE;	
}	

void xy_model::run(int N,int nsweeps,double J){
	/*
	 * Run function for the MCMC of xy model
	 * 
	 * Args
	 * N(int): The number of elements in each direction
	 * nsweeps(int): Number of sweeps to be performed 
	 * J(double): The energy constant in -J*cos(thetai-thetaj)
	 */

	omp_set_num_threads(8);
	// Initialized random lattice
	std::vector<std::vector<double>> lattice=Random_initialize(N);

	// Calculate initial energy
	double E = E_initial(lattice,J);

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			double updated_spin = rotatespin(lattice[i][j]); 
			double dE = calc_deltaE(lattice,updated_spin,i,j,J);
			std::cout << dE << std::endl;
		}
	}
}
