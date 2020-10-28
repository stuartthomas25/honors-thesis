#include <iostream>
#include <time.h>
#include "sweep.cpp"
#include <chrono>
#include <thread>
#include "mpi.h"
#include <algorithm>
#include <typeinfo>
#include <fstream>
#include <unistd.h>
#include <chrono>

using namespace std;
using namespace croutines;
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
    float m02;
    char* filename;
    while ((c = getopt(argc, (char **)argv, "o:m:")) != -1) {
        switch((char)c) {
            case 'o':    
                filename = optarg;
                break;
            case 'm':    
                m02 = atof(optarg);
                break;
            default:
                cout << "Unrecognized argument: " << (char)c << endl;
                return 1;

        }
    };
    int process_Rank, size_Of_Cluster;

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    if (DIM * DIM % size_Of_Cluster != 0) {
        cout << "ERROR: Size of lattice must be divisible by cluster size" << endl;
        return 1;
    }


    int* mpi_assignments;



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

    Sweeper sweeper(m02, 0.5, MPI_COMM_WORLD);

    auto start = high_resolution_clock::now();
    int sweeps = 10000;
    int thermalization = 1000;
    int record_rate = 100;
    vector<measurement> meas_data = sweeper.full_sweep(sweeps, thermalization, record_rate);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start); 

    avg_mag = 0;
    for (measurement m : meas_data) {
        outputfile << m.action << ", " << m.phibar[0] << endl;

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
