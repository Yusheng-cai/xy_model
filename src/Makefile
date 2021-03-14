# -*- MakeFile -*-
#
#  target: dependencies
#  	action

CC=/usr/local/bin/g++
OUTPUT=MCMC
INSTALLDIR=/usr/local/bin/

all:$(OUTPUT)

MCMC:main.o xy_model.o 
	$(CC) main.o xy_model.o -o MCMC -L/usr/local/lib -fopenmp 

main.o:main.cpp xy_model.h
	$(CC) -c main.cpp 

xy_model.o:xy_model.cpp
	$(CC) -c xy_model.cpp -fopenmp  

clean:
	rm -f *.o
	rm -f $(OUTPUT) 
