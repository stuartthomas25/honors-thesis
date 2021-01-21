export OMPI_CXX=g++-10
#export LDFLAGS="-L/usr/local/opt/binutils/lib"
export CPPFLAGS=-I/usr/local/opt/binutils/include
binary=bin/sweep
cppfiles=sweep.cpp phi.cpp main.cpp 
hfiles=sweep.h phi.h

.PHONY: clean
.DEFAULT_GOAL: sweep

sweep: $(cppfiles) $(hfiles)
	mpic++ $(cppfiles) -o $(binary) 

debug: main.cpp sweep.cpp sweep.h
	mpic++ -g $(cppfiles) -o $(binary) 

clean: 
	rm -f $(binary)
