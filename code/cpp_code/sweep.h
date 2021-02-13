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

#include "constants.h"
#include "mpi.h"
#include "phi.h"
#include "lattice.h"

using namespace std;

enum ClusterAlgorithm { WOLFF, SWENDSEN_WANG };

struct measurement {
    //Phi phibar;
    double action;
    double chi_m;
};

struct sweep_args {
    int sweeps;
    int thermalization;
    int record_rate;
    ClusterAlgorithm cluster_algorithm = WOLFF; 
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
    
    unordered_map<COLOR, vector<int>> mpi_assignments;

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
        vector<int> full_neighbors(int site);

    public:
        const int DIM;
        Lattice2D lat;
        double action;


        static double beta;
        Sweeper();
        ~Sweeper();
        Sweeper(int DIM, MPI_Comm c, bool makeGif);

        double full_action();
        void assert_action(double tol=0.001);
        int wrap(int c);
        static double lagrangian(Phi phi, Phi nphi_sum);
        double rand_dist(double r);
        Phi new_value(Phi old_phi);
        Phi proj_vec();
        Phi random_phi();

        vector<measurement> full_sweep(const sweep_args& args);
        tuple<vector<Phi>, double> sweep(COLOR color);
        void wolff();

        void broadcast_lattice();
        void collect_changes(vector<Phi> dphis, double dS, COLOR color);

        double Padd(Phi dphi, Phi phi_b);
        unordered_set<int> generate_cluster(int seed, Phi r, bool accept_all);

        void flow(double t); 

};


double randf();
int randint(int n);
Phi sign(Phi x);
double sign(double x);




#endif
