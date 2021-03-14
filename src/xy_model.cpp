#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "xy_model.h"
#include <map>
#include <omp.h>
#include <fstream>
#include <string>


#define PI 3.14159265
#define RAD 3.14159265/180.0
#define R 0.008314

double xy_model::randnum(double low, double high){
	/*
	 * Function that generates a uniformly generated random number between low and high
	 *
	 * Args:
	 * low(float): The lower bound of the random number 
	 * high(float): THe upper bound of the random number 
	 *
	 * Returns:
	 * (double): A random number between 0 and 1
	 */
	return low + (double) (high-low)*rand()/RAND_MAX;	
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
	double RANDANGLE = randnum(-180,180);
	
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
			lattice[i][j] = randnum(-180.0,180.0);
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

		std::cout << -J*(cos((currval-rval)*RAD)+cos((currval-lval)*RAD)+cos((currval-uval)*RAD)+\
				cos((currval-dval)*RAD)) << std::endl;	
		energy +=-J*(cos((currval-rval)*RAD)+cos((currval-lval)*RAD)+cos((currval-uval)*RAD)+\
				cos((currval-dval)*RAD));	
	}
	double end=omp_get_wtime();
	std::cout << end-start<<std::endl;

	return energy/2.0;
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
	
	double currspin = lattice[i][j];

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

	double dE =J*(cos((currspin-rval)*RAD) - cos((updated_spin-rval)*RAD)+\
			cos((currspin-lval)*RAD) - cos((updated_spin-lval)*RAD)+\
			cos((currspin-uval)*RAD) - cos((updated_spin-uval)*RAD)+cos((currspin-dval)*RAD) - cos((updated_spin - dval)*RAD));	
	
       return dE;	
}	

void xy_model::run(int N,int nsweeps,double J,double T,std::string ENERGY, std::string CONFIG){
	/*
	 * Run function for the MCMC of xy model
	 * 
	 * Args
	 * N(int): The number of elements in each direction
	 * nsweeps(int): Number of sweeps to be performed 
	 * J(double): The energy constant in -J*cos(thetai-thetaj)
	 * T(double): The temperature of the system
	 */
	
	omp_set_num_threads(8);
	// Initialized random lattice
	std::vector<std::vector<double>> lattice=Random_initialize(N);

	// create the files to be written to 
	std::ofstream energyfile;
	std::ofstream configfile;
	energyfile.open(ENERGY,std::ios::out);
	configfile.open(CONFIG,std::ios::out);

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			configfile << lattice[i][j] << "\t";
		}
		configfile << "\n";
	}


	// Calculate initial energy
	double E = E_initial(lattice,J);
	std::cout << "Initial E is " << E << std::endl;
	double RT = R*T;	
	
	for(int m=0;m<nsweeps;m++){
		configfile << "nsweep  " << m+1 << "\n";
		energyfile << "nsweep  " << m+1 << "\n";
		for(int k=0;k<N*N;k++){
			int num = (int)round(randnum(0.0,1.0)*(N*N-1));
			int i=num/N;
			int j=num%N;
			
			double updated_spin = rotatespin(lattice[i][j]); 
			double dE = calc_deltaE(lattice,updated_spin,i,j,J);
			
			// if dE is less than 0, then accept the move
			if (dE < 0){
				lattice[i][j] = updated_spin;
				configfile << i << "\t";
				configfile << j << "\t";
				configfile << lattice[i][j] << "\n";

				E = E + dE;
				// std::cout << "Accepted because dE is negative and dE: " << dE << std::endl;
			}
			else{
				double factor=exp(-dE/(RT));
				double r = randnum(0.0,1.0);
				
				if (r < factor){
					lattice[i][j] = updated_spin;
					configfile << i << "\t";
					configfile << j << "\t";
					configfile << lattice[i][j] << "\n";


					E = E + dE;
					// std::cout << "dE is positive, dE is " << dE << "factor is " << factor << "randnum is " << r << std::endl;
				}
				else{
					configfile << i << "\t";
					configfile << j << "\t";
					configfile << "Unchanged" << "\n";	
				}
			}
			energyfile << E;
			energyfile << "\n";
		}
	}
	energyfile.close();
	configfile.close();
	std::cout << "Final energy is " << E << std::endl;
}
