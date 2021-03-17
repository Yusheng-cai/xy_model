# -*- MakeFile -*-
#
#  target: dependencies
#  	action

GCC=/usr/local/bin/g++
OUTPUT=xy
INSTALLDIR=/usr/local/bin/

all:$(OUTPUT)

MCMC:main.o xy_model.o 
	$(GCC) -std=c++11 main.o xy_model.o -o $(OUTPUT) -L/usr/local/lib -fopenmp 

main.o:main.cpp xy_model.h
	$(GCC) -std=c++11 -c main.cpp 

xy_model.o:xy_model.cpp
	$(GCC) -std=c++11 -c xy_model.cpp -fopenmp  

clean:
	rm -f *.o
	rm -f $(OUTPUT) 

install:
	cp $(OUTPUT) $(INSTALLDIR)
