export OMPI_CXX=g++-10
#export LDFLAGS="-L/usr/local/opt/binutils/lib"
export CPPFLAGS=-I/usr/local/opt/binutils/include
binaries=bin/sweep

.PHONY: clean
.DEFAULT_GOAL: sweep

sweep: main.cpp sweep.cpp sweep.h
	mpic++ -o bin/sweep main.cpp

profile: main.cpp sweep.cpp sweep.h
	mpic++ $(CPPFLAGS) -o bin/sweep -pg main.cpp

#debug: main.cpp sweep.cpp sweep.h
	#mpic++ -g -o bin/sweep main.cpp

clean: 
	rm -f $(binaries)
