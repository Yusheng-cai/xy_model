# -*- MakeFile -*-
#
#  target: dependencies
#  	action

CC=/usr/local/bin/g++

all:MCMC

MCMC:main.o simulation.o 
	$(CC) main.o simulation.o -o MCMC -L/usr/local/lib 

main.o:main.cpp simulation.h
	$(CC) -c main.cpp 

simulation.o:simulation.cpp
	$(CC) -c simulation.cpp
