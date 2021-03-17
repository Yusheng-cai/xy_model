#include "xy_model.h"

#define PI 3.14159265
#define RAD 3.14159265/180.0

double xy_model::randnum(double low, double high){
	/*
	 * Function that generates a uniformly generated random number between low and high
	 *
	 * Args:
	 * low(double): The lower bound of the random number 
	 * high(double): THe upper bound of the random number 
	 *
	 * Returns:
	 * (double): A random number between 0 and 1
	 */
	return low + (double) (high-low)*rand()/RAND_MAX;	
}

double xy_model::rotatespin(const double& angle){
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

double xy_model::E_initial(const std::vector<std::vector<double>>& lattice){
	/*
	 * Function that calculates the energy of the initial spin system E/J
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

		energy += -(cos((currval-rval)*RAD)+cos((currval-lval)*RAD)+cos((currval-uval)*RAD)+\
				cos((currval-dval)*RAD));	
	}

	return energy/2.0;
}

double xy_model::calc_deltaE(const std::vector<std::vector<double>> &lattice, double updated_spin,int i, int j){
	/*
	 * Function that calculates the unitless change of energy from changing the spin dE/J
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
	int rows=lattice.size();
	int cols=lattice[0].size(); 
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

	double dE =cos((currspin-rval)*RAD) - cos((updated_spin-rval)*RAD)+\
			cos((currspin-lval)*RAD) - cos((updated_spin-lval)*RAD)+\
			cos((currspin-uval)*RAD) - cos((updated_spin-uval)*RAD)+\
			cos((currspin-dval)*RAD) - cos((updated_spin - dval)*RAD);	
	
       return dE;	
}	

void xy_model::run(int N,int nsweeps,double kbTJ,double r,std::string EDR, std::string XTC,std::string P,int printevery){
	/*
	 * Run function for the MCMC of xy model
	 * 
	 * Args
	 * N(int): The number of elements in each direction
	 * nsweeps(int): Number of sweeps to be performed 
	 * kbTJ(double): kbT/J
	 * r(double): The raio of sweeps to start performing averaging
	 * EDR(std::string): A string for the name of the energy file output
	 * XTC(std::string): A string for the name of the configuration file output
	 * P(std::string): A string for the name of the property file output
	 * printevery(int): The frequency of iterations at which to write data to energy file
	 */
	
	omp_set_num_threads(8);

	// updated set to false initially
	bool updated = 0;
	double iteration = 1.0;
	double cosaverage = 0.0;
	double sinaverage = 0.0;
	double N2 = N*N;


	// Initialized random lattice 
	std::vector<std::vector<double>> lattice=Random_initialize(N);

	// create the files to be written to 
	std::ofstream energyfile;
	std::ofstream configfile;
	std::ofstream propertyfile;
	energyfile.open(EDR,std::ios::out);
	configfile.open(XTC,std::ios::out);	
	propertyfile.open(P,std::ios::out);
	energyfile << "Energy(unitless)\tAverage Cosine\tAverage Sine\n";


	// Calculate initial energy unitless
	double E = E_initial(lattice);
	double E_avg = E;
	double E2_avg = E*E;

	std::cout << "Initial unitless energy E/J is " << E << std::endl;

	// Calculate the average cosine and sin theta	
	#pragma omp parallel for shared(lattice) reduction(+:cosaverage,sinaverage)	
	for(int k=0;k<N*N;k++){
		int i=k/N;
		int j=k%N;
		double c = cos(lattice[i][j]*RAD)/N2;
		double s = sin(lattice[i][j]*RAD)/N2;
		cosaverage += c; 
		sinaverage += s;
	}

	
	for(int m=0;m<nsweeps;m++){
		for(int k=0;k<N*N;k++){
			updated = 0;
			// old lattice value
			double old=0;

			// random number to decide which lattice to perform calculation on
			int num = (int)round(randnum(0.0,1.0)*(N*N-1));
			int i=num/N;
			int j=num%N;
		
			// update spin and calculate change in energy	
			double updated_spin = rotatespin(lattice[i][j]); 
			double dE = calc_deltaE(lattice,updated_spin,i,j);
			
			// if dE is less than 0, then accept the move, else accept with MC criteria 
			if (dE < 0){
				old = lattice[i][j];
				lattice[i][j] = updated_spin;
				E = E + dE;

				updated = 1;
			}
			else{
				double factor=exp(-dE/(kbTJ));
				double r = randnum(0.0,1.0);
				
				if (r < factor){
					old = lattice[i][j];
					lattice[i][j] = updated_spin;
					E = E + dE;

					updated = 1;
				}
				
			}
			// Update average energy and variance
			if (m >= (double)nsweeps*r){
				E2_avg  = (E2_avg*iteration + E*E)/(iteration+1); 
				E_avg = (E_avg*iteration + E)/(iteration + 1);
				iteration += 1;
			}

				
			// write to energy file every printevery
			if((int)iteration %printevery == 0){
				if(updated == 1){
					double cosold = cos(old*RAD);
					double cosnew = cos(lattice[i][j]*RAD);
					double sinold = sin(old*RAD);
					double sinnew = sin(lattice[i][j]*RAD);

					cosaverage = (cosaverage*N2 - cosold + cosnew)/N2;
					sinaverage = (sinaverage*N2 - sinold + sinnew)/N2;
				}
				energyfile << iteration << "\t" << E <<"\t" << cosaverage << "\t" << sinaverage << "\n";
			}
		}

		// store coordinates every 10 sweeps
		if(m % 10 == 0){
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					configfile << lattice[i][j] << "\t";
				}
				configfile << "\n";
			}
			}
}

	double variance = E2_avg - pow(E_avg,2);
	energyfile.close();
	configfile.close();

	propertyfile<< "Average energy (unitless): " << E_avg << "\n";
	propertyfile<< "Variance : " << variance << "\n";
	propertyfile<< "Heat capacity (kbT): " << variance/(N*N*pow(kbTJ,2)) << "\n";
	propertyfile.close();
}
