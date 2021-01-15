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
#include <map>
#include <string>
#include <memory>
#include "mpi.h"
#include "gif.h"

using namespace std;

#define N 1
#define MASTER 0
#define PROG_CHAR "#"

namespace croutines {

    enum ClusterAlgorithm { WOLFF, SWENDSEN_WANG };
    
    class Phi {
        array<double, N> phi;

        public:
            Phi();
            Phi(array<double, N> phi);
            void init_as_zero();
            double norm_sq() const;
            Phi operator+ (const Phi & phi) const;
            Phi operator- (const Phi & phi) const;
            Phi operator- () const;
            double operator* (const Phi & phi) const;
            Phi operator* (const double & a) const;
            friend ostream& operator<< (ostream& os, const Phi & aPhi);
            double operator[] (int i) const;
            double& operator[] (int i);
            bool operator== (const Phi & phi) const;
            void print();
    };
      


    struct measurement {
        Phi phibar;
        double action;
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
        const int DIM;

        enum COLOR {
            black,
            white
        };
        
        map<COLOR, vector<int>> mpi_assignments;

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
            vector<uint8_t> gif_vec;
            GifWriter gif_writer;
            auto static constexpr gif_filename = "lattice.gif";
            int gif_delay = 10;

        public:
            int dim;
            vector<Phi> lat;
            double action;


            double redef_mass, quarter_lam, beta;
            Sweeper();
            ~Sweeper();
            Sweeper(double m02, double lam, int DIM, MPI_Comm c, bool makeGif);

            void write_gif_frame(double rate);

            double full_action();
            int wrap(int c);
            void full_neighbors(int site, int neighbors[4]);
            double lagrangian(Phi phi, Phi nphi_sum);
            double phi4_lagrangian(Phi phi, Phi nphi_sum);
            double sigma_lagrangian(Phi phi, Phi nphi_sum);
            double rand_dist(double r);
            Phi new_value(Phi old_phi);
            Phi proj_vec();
            Phi random_phi();

            vector<measurement> full_sweep(const sweep_args& args);
            tuple<vector<Phi>, double> sweep(COLOR color);
            void wolff();

            void broadcast_lattice();
            void collect_changes(vector<Phi> dphis, double dS, COLOR color);

            double Padd(Phi phi_a, Phi phi_b);
            set<int> generate_cluster(int seed, bool accept_all);

    };


    double randf();
    int randint(int n);
    Phi sign(Phi x);
    double sign(double x);




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
            //os << '\n';
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
