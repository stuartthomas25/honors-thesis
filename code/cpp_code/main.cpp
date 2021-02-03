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
#include "lattice.cpp"

using namespace std;
using namespace std::chrono;




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





   /* vector<Phi> data;*/
    //for (int i=0; i<pow(DIM, 2); i++) {
        //if (process_Rank == 0) {
            //data.push_back( Phi({ rand_dist(randf()) }) );
        //} else {
            //data.push_back(Phi({0}));
        //}
    /*}*/

    double action;

    srand(time(NULL)); rand();

    vector<Phi> output;
    double dS, avg_mag;

    //tie(output, dS) = sweeper.sweep(data, sites);
    
    ofstream outputfile;
    outputfile.open(filename);


    Sweeper sweeper(dim, MPI_COMM_WORLD, gif);

    auto start = high_resolution_clock::now();
    int measurements = 100;
    int thermalization = 1000;
    int record_rate = 100;
    int sweeps = record_rate * measurements + thermalization;

    sweep_args args {
        .sweeps = sweeps,
        .thermalization = thermalization,
        .record_rate = record_rate,
        .progress = progress
    };
    vector<measurement> meas_data = sweeper.full_sweep(args);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start); 
    
    if (args.progress) {
        cout << meas_data.size() << " measurements" << endl;
    }
    avg_mag = 0;
    for (measurement m : meas_data) {
        outputfile << m.action;
        for (int i=0; i<N; i++) {
            outputfile << ", " << m.phibar[i];
        }
        outputfile << endl;
    }

    outputfile.close();

    cout << '\r' << "Wrote to " << filename << " in " << duration.count() << "s   \n\n";


   /* for (int i=0; i<4; i++) {*/
        /*for (int j=0; j<4; j++) {*/
            /*cout << sweeper.lat[i+DIM*j] << "  ";*/
        /*}*/
        /*cout << endl;*/
    /*}*/
    //cout << meas_data.size() << endl;
    //for (measurement md : meas_data) {
        //cout << md.action << " ";
    //}
    //cout << endl;

    //cout << dS << endl;
    //for(int i=0; i<dim; i++) {
        //for(int j=0; j<dim; j++) {
            //Phi phi = output[i*dim + j];
            //if (phi[0] < 0) {
                //cout << "-";
            //}else{
                //cout << "*";
            //}
        //}
        //cout << endl;
    MPI_Finalize();
    return 0;
};