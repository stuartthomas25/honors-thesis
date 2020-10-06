#ifndef SWEEP_H
#define SWEEP_H
#include <vector>
#include <tuple>
#include <stack>
#include <set>
#include <limits>
#include <array>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <ostream>
#include <string>
#include "mpi.h"

using namespace std;

#define N 2
#define PROG_CHAR "#"

namespace croutines {

    enum ClusterAlgorithm { WOLFF, SWENDSEN_WANG };

    class Phi {
        array<double, N> phi;

        public:
            Phi();
            Phi(array<double, N> phi);
            void init_as_zero();
            Phi operator+ (const Phi & phi) const;
            Phi operator- (const Phi & phi) const;
            Phi operator- () const;
            double operator* (const Phi & phi) const;
            Phi operator* (const double & a) const;
            double operator[] (int i) const;
            double& operator[] (int i);
            bool operator== (const Phi & phi) const;
            void print();
    };
            
    struct measurement {
        Phi phibar;
        double action;
    };

    class Sweeper {
        public:
            int dim;
            vector<Phi> lat;
            double redef_mass, quarter_lam, beta;
            vector<Phi> lattice;
            Sweeper();
            Sweeper(vector<Phi> lattice, double m02, double lam);
            int wrap(int c);
            void full_neighbors(int site, int neighbors[4]);
            double lagrangian(Phi phi, Phi nphi_sum);
            double phi4_lagrangian(Phi phi, Phi nphi_sum);
            double sigma_lagrangian(Phi phi, Phi nphi_sum);
            double rand_dist(double r);
            Phi new_value(Phi old_phi);
            Phi proj_vec();

            vector<measurement> full_sweep(int sweeps, int thermalization, int record_rate, ClusterAlgorithm cluster_algorithm, int cluster_rate);
            tuple<vector<Phi*>, double> sweep(vector<int>& sites);
            tuple<vector<Phi>, double> wolff();

            double Padd(Phi phi_a, Phi phi_b);
            set<int> generate_cluster(int seed, bool accept_all);

    };


    double randf();
    int randint(int n);
    Phi sign(Phi x);




    class progress_bar
    {
        static const auto overhead = sizeof " [100%]";

        ostream& os;
        const size_t bar_width;
        string message;
        const string full_bar;

     public:
        progress_bar(ostream& os, size_t line_width,
                     string message_, const char symbol = '.')
            : os{os},
              bar_width{line_width - overhead},
              message{move(message_)},
              full_bar{string(bar_width, symbol) + string(bar_width, ' ')}
        {
            if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
                os << message << '\n';
                message.clear();
            } else {
                message += ' ';
            }
            write(0.0);
        }

        // not copyable
        progress_bar(const progress_bar&) = delete;
        progress_bar& operator=(const progress_bar&) = delete;

        ~progress_bar()
        {
            write(1.0);
            os << '\n';
        }

        void write(double fraction);
    };

    void progress_bar::write(double fraction)
    {
        // clamp fraction to valid range [0,1]
        if (fraction < 0)
            fraction = 0;
        else if (fraction > 1)
            fraction = 1;

        auto width = bar_width - message.size();
        auto offset = bar_width - static_cast<unsigned>(width * fraction);

        os << '\r' << message;
        os.write(full_bar.data() + offset, width);
        os << " [" << setw(3) << static_cast<int>(100*fraction) << "%] " << flush;
    } 

}

#endif
