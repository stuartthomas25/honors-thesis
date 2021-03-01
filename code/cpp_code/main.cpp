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
#include "Yaml.hpp"
#include "sweep.h"

#include "observables.h"

using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[]) {

    int c;
    char* filename;
    char* in_filename;
    bool progress = true;

    while ((c = getopt(argc, (char **)argv, "qi:o:")) != -1) {
        switch((char)c) {
            case 'i':    
                in_filename = optarg;
                break;
            case 'o':    
                filename = optarg;
                break;
            case 'q':
                progress = false;
                break;
            default:
                cout << "Unrecognized argument: " << (char)c << endl;
                return 1;

        }
    };

    Yaml::Node root;
    Yaml::Parse(root, in_filename);
    
    double beta = root["beta"].As<double>();
    int dim = root["L"].As<int>();
    int measurements = root["measurements"].As<int>();
    int thermalization = root["thermalization"].As<int>();
    int record_rate = root["record_rate"].As<int>();
    int cluster_rate = 5;

    vector<double> taus;
    auto iter =  root["taus"].Begin();
    while (iter != root["taus"].End()) {
        taus.push_back((*iter).second.As<double>());
        iter++;
    }


    Lattice2D::L = dim;
    Lattice2D::beta = beta;
    int process_Rank, size_Of_Cluster;

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    if (dim * dim % size_Of_Cluster != 0) {
        cout << "ERROR: Size of lattice must be divisible by cluster size" << endl;
        return 1;
    }

    srand(time(NULL)); rand();

    vector<Phi> output;
    double dS, avg_mag;

    //tie(output, dS) = sweeper.sweep(data, sites);
    
    Sweeper sweeper(dim, MPI_COMM_WORLD);

    auto start = high_resolution_clock::now();
    int sweeps = record_rate * measurements + thermalization;

    vector<BaseObservable*> observables{
        new observables::beta(),
        new observables::L(),
        //new observables::chi_m(),
        new observables::action(), 
        new observables::Q()
    };

    Recorder recorder(observables);
    recorder.reserve(measurements);

    const sweep_args args = {
        sweeps,
        thermalization,
        record_rate,
        ClusterAlgorithm::WOLFF,
        taus,
        cluster_rate,
        progress
    };

    sweeper.full_sweep(&recorder, args);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start); 
    
    if (args.progress) {
        cout << recorder.size() << " measurements" << endl;
    }
    recorder.write(filename);

    cout << '\r' << "Wrote to " << filename << " in " << duration.count() << "s   \n\n";
    MPI_Finalize();
    return 0;
};
