export OMPI_CXX=g++

test: main.cpp sweep.cpp sweep.h
	mpic++ -o bin/sweep main.cpp

debug: main.cpp sweep.cpp sweep.h
	mpic++ -g -o bin/sweep main.cpp
