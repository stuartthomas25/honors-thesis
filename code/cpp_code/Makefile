export OMPI_CXX=g++
binaries=bin/sweep

test: main.cpp sweep.cpp sweep.h
	mpic++ -o bin/sweep main.cpp

#debug: main.cpp sweep.cpp sweep.h
	#mpic++ -g -o bin/sweep main.cpp

#clean: 
	#rm -f $(binaries)
