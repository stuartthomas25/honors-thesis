#include <iostream>
#include <time.h>
#include <chrono>
#include <thread>
#include "mpi.h"
#include <algorithm>
#include <typeinfo>
#include <fstream>
#include <unistd.h>
#include <chrono>

#include "phi.h"
#include "sweep.h"

using namespace std;
using namespace std::chrono;


/*void update_lattice(vector<Phi*>* dphis, double dS, int master = 0) {*/
    //int process_Rank, size_Of_Cluster;
    //MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    //MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    //const int head = 0 ;

    //const int raw_data_len = DIM*DIM*N;
    //double raw_data[raw_data_len];
    //if (process_Rank == master) {
        //for (int i = 0; i < raw_data_len; i++) {
            //raw_data[i] = (*data)[i/N] [i%N];
        //};
    //}
    



/*}*/



int main(int argc, char *argv[]) {

    int c;
    float beta;
    char* filename;
    int dim;
    bool progress = true;
    bool gif = false;
    while ((c = getopt(argc, (char **)argv, "gqo:b:L:")) != -1) {
        switch((char)c) {
            case 'o':    
                filename = optarg;
                break;
            case 'b':    
                beta = atof(optarg);
                break;
            case 'L':    
                dim = atof(optarg);
                break;
            case 'q':
                progress = false;
                break;
            case 'g':
                gif = true;
                break;
            default:
                cout << "Unrecognized argument: " << (char)c << endl;
                return 1;

        }
    };
    Lattice2D::L = dim;
    Sweeper::beta = beta;
    int process_Rank, size_Of_Cluster;

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    if (dim * dim % size_Of_Cluster != 0) {
        cout << "ERROR: Size of lattice must be divisible by cluster size" << endl;
        return 1;
    }


    double action;

    srand(time(NULL)); rand();

    vector<Phi> output;
    double dS, avg_mag;

    //tie(output, dS) = sweeper.sweep(data, sites);
    
    ofstream outputfile;
    outputfile.open(filename);

    Sweeper sweeper(dim, MPI_COMM_WORLD, gif);

    sweeper.flow(10);

    outputfile.close();
    MPI_Finalize();
    return 0;
};
