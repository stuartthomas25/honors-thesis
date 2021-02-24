export OMPI_CXX=/usr/bin/clang++
#export LDFLAGS="-L/usr/local/opt/binutils/lib"
export CPPFLAGS=-I/usr/local/opt/binutils/include
binary=bin/sweep
classes=sweep.cpp phi.cpp lattice.cpp Yaml.cpp
hfiles=sweep.h phi.h
ifdef debug
	options=-std=c++17 -Ofast -g -fprofile-instr-generate
else 
	options=-std=c++17 -Ofast
endif 

export PATH="/opt/homebrew/bin:$PATH"
mpi_path=/opt/homebrew/bin/mpic++

.PHONY: clean
.DEFAULT_GOAL: sweep

sweep: main.cpp $(classes) $(hfiles)
	$(mpi_path) main.cpp $(options) $(classes) -o $(binary) 

asm: main.cpp $(classes) $(hfiles)
	$(mpi_path)  main.cpp $(options) -S $(classes)
	wc -l *.s

flow: flow.cpp $(classes) $(hfiles)
	$(mpi_path)  $(options) flow.cpp $(classes) -o bin/flow


clean: 
	rm -f $(binary)
