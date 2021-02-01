export OMPI_CXX=g++-10
#export LDFLAGS="-L/usr/local/opt/binutils/lib"
export CPPFLAGS=-I/usr/local/opt/binutils/include
binary=bin/sweep
classes=sweep.cpp phi.cpp lattice.cpp
hfiles=sweep.h phi.h
options=-std=c++20 -g


.PHONY: clean
.DEFAULT_GOAL: sweep

sweep: main.cpp $(classes) $(hfiles)
	mpic++ main.cpp $(options) $(classes) -o $(binary) 

flow: flow.cpp $(classes) $(hfiles)
	mpic++ $(options) flow.cpp $(classes) -o bin/flow 


clean: 
	rm -f $(binary)
