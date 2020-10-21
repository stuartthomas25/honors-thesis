#include <iostream>
#include <time.h>
#include "sweep.cpp"
#include <chrono>
#include <thread>
#include "mpi.h"
#include <algorithm>
#include <typeinfo>
#include <fstream>


using namespace std;
using namespace croutines;


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



int main() {

    int process_Rank, size_Of_Cluster;

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    if (DIM * DIM % size_Of_Cluster != 0) {
        cout << "ERROR: Size of lattice must be divisible by cluster size" << endl;
        return 1;
    }


    int* mpi_assignments;


    vector<double> m02s = {-0.80, -0.76, -0.72, -0.68, -0.64};

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
    outputfile.open("actions.csv");
    for (double m02 : m02s) {
        Sweeper sweeper(m02, 0.5, MPI_COMM_WORLD);

        vector<measurement> meas_data = sweeper.full_sweep(10000, 1000, 10);
        
        avg_mag = 0;
        for (measurement m : meas_data) {
            avg_mag += abs(m.phibar[0]);
        }
        avg_mag /= meas_data.size();

        outputfile << 0 << ", " << avg_mag << endl;
    }
    outputfile.close();


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
