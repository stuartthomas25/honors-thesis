#ifndef SWEEP_H
#define SWEEP_H
#include <vector>
#include <tuple>
#include <stack>
#include <unordered_set>
#include <limits>
#include <array>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <ostream>
#include <memory>
#include <string>

#include "constants.h"
#include "mpi.h"
#include "phi.h"
#include "lattice.h"


using namespace std;

enum ClusterAlgorithm { WOLFF, SWENDSEN_WANG };
typedef int site;


class Recorder {
    typedef double (*ObservableFunc) (Lattice2D& lat);
    vector<double> measurements;
    vector<ObservableFunc> observables;
    public:
        Recorder(vector<ObservableFunc> some_observables);
        void reserve(size_t size);
        int size();
        void record(Lattice2D& lat, double x);
        void write(string filename);
};


struct sweep_args {
    int sweeps;
    int thermalization;
    int record_rate;
    ClusterAlgorithm cluster_algorithm = WOLFF; 
    vector<double> ts;
    int cluster_rate = 5;
    bool progress = true;
};

class Sweeper {
    const int process_Rank;
    const int size_Of_Cluster;
    const int sites_per_node;

    enum COLOR {
        black,
        white
    };
    
    unordered_map<COLOR, vector<site>*> mpi_assignments;

    private:
        static int get_rank(MPI_Comm c) {
            int rank;
            MPI_Comm_rank(c, &rank);
            return rank;
        }
        static int get_size(MPI_Comm c) {
            int size;
            MPI_Comm_size(c, &size);
            return size;
        }

        bool gif;
        auto static constexpr gif_filename = "lattice.gif";
        int gif_delay = 10;
        vector<site> full_neighbors(site aSite);
        Plaquette plaquette(site aSite);
        vector<Phi> dphis;
        Lattice2D flowed_lat;
        Lattice2D k1, k2, k3, k4;

    public:
        const int DIM;
        Lattice2D lat;


        static double beta;
        Sweeper();
        ~Sweeper();
        Sweeper(int DIM, MPI_Comm c);

        double full_action();
        void assert_action(double tol=0.001);
        int wrap(site c);
        static double lagrangian(Phi phi, Phi nphi_sum);
        double rand_dist(double r);
        Phi new_value(Phi old_phi);
        Phi proj_vec();
        Phi random_phi();

        void full_sweep(Recorder* recorder, const sweep_args& args);
        double sweep(COLOR color);
        void wolff();

        void broadcast_lattice();
        void collect_changes(double dS, COLOR color);

        double Padd(Phi dphi, Phi phi_b);
        unordered_set<int> generate_cluster(int seed, Phi r, bool accept_all);

        void flow(vector<double> ts, Recorder* recorder=nullptr); 

};


double randf();
int randint(int n);
Phi sign(Phi x);
double sign(double x);




#endif
