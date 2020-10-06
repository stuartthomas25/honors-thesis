#include <iostream>
#include <time.h>
#include "sweep.cpp"
#include <chrono>
#include <thread>
#include "mpi.h"
#include <algorithm>
#include <typeinfo>

#define DIM 4

using namespace std;
using namespace croutines;

void broadcast_lattice(vector<Phi>* data, int master=0) {
    int process_Rank, size_Of_Cluster;
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    int head = 0;

    const int raw_data_len = DIM*DIM*N;
    double raw_data[raw_data_len];
    if (process_Rank == master) {
        for (int i = 0; i < raw_data_len; i++) {
            raw_data[i] = (*data)[i/N] [i%N];
        };
    }


    MPI_Bcast(&raw_data, raw_data_len, MPI_DOUBLE, master, MPI_COMM_WORLD);


    if (process_Rank != master) {
        for (int i = 0; i < raw_data_len; i++) {
            (*data)[i/N] [i%N] = raw_data[i];
        };
    }






}

int main() {
    int process_Rank, size_Of_Cluster;

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);


    vector<Phi> data;
    for (int i=0; i<pow(DIM, 2); i++) {
        if (process_Rank == 0) {
            data.push_back(Phi({-(double)i,(double)i}));
        } else {
            data.push_back(Phi({0,0}));
        }
    }

    cout << process_Rank << " 1 " << data[5][1] << endl;
    broadcast_lattice(&data);
    cout << process_Rank << " 2 " << data[5][1] << endl;

    MPI_Finalize();
    return 0;

    //srand(time(NULL)); rand();

    //vector<Phi> data;
    //for (int i=0; i<pow(DIM, 2); i++) {
        //data.push_back(Phi({-1}));
    //}
    //vector<int> sites = {0, 1, 2, 3};

    //vector<Phi> output;
    //double dS;

    ////tie(output, dS) = sweeper.sweep(data, sites);
    

    /*Sweeper sweeper(data, -0.68, 0.5);*/

    /*vector<measurement> meas_data = sweeper.full_sweep(1000, 10, 1);*/
    //cout << meas_data.size() << endl;
    //for (measurement md : meas_data) {
        //cout << md.action << " ";
    //}
    //cout << endl;

    //[>cout << dS << endl;<]
    ////for(int i=0; i<dim; i++) {
        ////for(int j=0; j<dim; j++) {
            ////Phi phi = output[i*dim + j];
            ////if (phi[0] < 0) {
                ////cout << "-";
            ////}else{
                ////cout << "*";
            ////}
        ////}
        ////cout << endl;
    //[>}<]

    /*return 0;*/
};
