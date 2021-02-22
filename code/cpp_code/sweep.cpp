#include <iostream>
#include <stdexcept>
#include "sweep.h"
#include <algorithm>
#include "gif.h"
#include "progress.cpp"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

//#define VERIFY_ACTION
//#define GIF
//#define ACTIONFLOW
#define ADAPTIVESTEP

using namespace std;

Recorder::Recorder(vector<BaseObservable*> some_observables) {
    observables = some_observables;
}

void Recorder::record(Lattice2D& lat, double x=0) {
    measurements.push_back( x );
    for (const BaseObservable* f : observables) {
        measurements.push_back( (*f)(lat) );
    };
}

void Recorder::reserve(size_t size) {
    return measurements.reserve((observables.size()+1) * size);
}

int Recorder::size() {
    return measurements.size() / (observables.size() + 1);
}

void Recorder::write(string filename) {
    ofstream outputfile;
    outputfile.open(filename);
    outputfile << "tau";
    for (const auto* obs : observables) {
        outputfile << "," << obs->name();
    }
    outputfile << endl;
    

    for (int i=0; i<measurements.size(); i++) {
        if ( (i+1)%(observables.size()+1) != 0 ) {
            outputfile << measurements[i] << ",";
        } else {
            outputfile << measurements[i] << endl;
        };
    }
    outputfile.close();
}

Recorder::~Recorder() {
    for (const BaseObservable* f : observables) {
        delete f;
    };

}


Sweeper::Sweeper (int DIM, MPI_Comm c) : 
    process_Rank(get_rank(c)), 
    size_Of_Cluster(get_size(c)), 
    sites_per_node((DIM*DIM)/2 / get_size(c)),
    DIM{DIM}
{
    // generate lattice
    array<double, N> zero_array = {};
    Phi zero_phi(zero_array);

    dphis.resize(sites_per_node);

    mpi_assignments.push_back(new vector<site>);
    mpi_assignments.push_back(new vector<site>);

    site s;
    for (int i = 0; i<DIM; i++){
        for (int j = 0; j<DIM; j++) {
            s = i*DIM + j;
            lat[s] = new_value(zero_phi);
            
            lat.neighbor_map.push_back(full_neighbors(s)); 
            lat.plaquette_map.push_back(plaquette(s)); 

            // generate MPI assignments
            if (s / (2*sites_per_node) == process_Rank) { // Factor of 2 for checkerboard
                if ( (i+j)%2 ) {
                    mpi_assignments[COLOR::black]->push_back(s);
                } else {
                    mpi_assignments[COLOR::white]->push_back(s);
                }   
            }
        }
    }

    lat.action = full_action(); 

    int offset = process_Rank *  sites_per_node;

}

Sweeper::~Sweeper() {
    delete mpi_assignments[COLOR::black];
    delete mpi_assignments[COLOR::white];

}

void print(std::vector <Phi> const &a) {
   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
   cout << endl;
}

double Sweeper::full_action() {
    site s;
    int i;
    Phi forward_nphi_sum;
    vector<site> neighbors;
    double S = 0;
    // initial action
    for (int s=0; s<DIM*DIM; s++) {
        neighbors = lat.neighbor_map[s];
        forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

        S += lat.lagrangian(lat[s], forward_nphi_sum);
    }
    return S;

}

int Sweeper::wrap( int c) {
    int mod = c % DIM;
    if (mod < 0) {
        mod += DIM;
    }
    return mod;

}

vector<site> Sweeper::full_neighbors(site aSite) {
    int x = aSite % DIM;
    int y = aSite / DIM;
    vector<site> neighbors; 
    neighbors.push_back(wrap(x-1) + DIM * y);
    neighbors.push_back(x + DIM * wrap(y-1));
    neighbors.push_back(wrap(x+1) + DIM * y);
    neighbors.push_back(x + DIM * wrap(y+1));
    return neighbors;
}

Plaquette Sweeper::plaquette(site aSite) {
    int x = aSite % DIM;
    int y = aSite / DIM;
    return make_tuple(  x + DIM * y,
                        wrap(x+1) + DIM * y,
                        wrap(x+1) + DIM * wrap(y+1),
                        x + DIM * wrap(y+1)
                     );
}


double Sweeper::rand_dist(double r) {
    return 3 * r - 1.5;
}

/*    Phi Sweeper::new_value(Phi old_phi) {*/


    //Phi dphi;
    //double sum_of_squares = 0;
    //for (int i = 0; i<N-1; i++) {
        //dphi[i] = rand_dist(randf());

        //sum_of_squares += pow(old_phi[i] + dphi[i], 2);
    //}
    //dphi[N-1] = (2*randint(2) - 1) * sqrt(1-sum_of_squares) - old_phi[N-1]; //the first term randomizes the sign

    //// test?
    //cout << "Dot: " << (old_phi + dphi) * (old_phi + dphi) << endl;
    //return dphi;
/*}*/

Phi Sweeper::new_value(Phi old_phi) {
    //Phi dphi;
    //dphi[0] = rand_dist(randf());
    //return dphi;
    return random_phi();
}

Phi Sweeper::proj_vec () {
    //for (int i = 0; i<N-1; i++) {
        //r[i] = 
    //}
    return Phi();
}

Phi Sweeper::random_phi() {
    Phi new_phi;
    for (int i=0; i<N; i++) {
        new_phi[i] = 2 * randf() - 1;
        if (abs(new_phi[i]) > 1e+10) {
            cout << "overflow" << endl;
            exit(1);

        }
    }
    return new_phi * (1/sqrt(new_phi.norm_sq()));

}

void write_gif_frame(Lattice2D& lat, GifWriter* gif_writer, int delay, double rate=3) {
    Phi aPhi;
    auto DIM = lat.L;
    vector<uint8_t> vec(4*DIM*DIM);

    for (int i=0; i<DIM*DIM; i++) {
        aPhi = lat[i];
        for (int j=0; j<N; j++) {
            vec[4*i + j] = static_cast<uint8_t>((aPhi[j]+1)*128);
        }
        vec[4*i + 3] = 255;
        //cout << vec[0] <<" "<< vec[1] <<" "<< vec[2] <<" "<< vec[3] <<" "<< endl;
    }
    GifWriteFrame(gif_writer, vec.data(), DIM, DIM, delay);
}

void Sweeper::assert_action(double tol) {
    double fa = full_action();
    if (abs(fa-lat.action)>tol) {
        cout << "ASSERT ACTION FAILED (" << lat.action << " != " << fa << ")\n";
        exit(1);
    } else {
        cout << "assert action passed (" << lat.action << " == " << fa << ")\n";

    }
}


void Sweeper::full_sweep(Recorder* recorder, const sweep_args& args = sweep_args()) {
   
    //if (args.cluster_algorithm != WOLFF) throw invalid_argument("Currently, Wolff is the only allowed algorithm");

    shared_ptr<progress_bar> progress;
    progress = args.progress ? make_shared<progress_bar>(cout, 70u, "Working") :  nullptr; 

    double dS;
    Phi phibar;
    double chi_m;

    int s;


    vector<uint8_t> white_vec(DIM*DIM*4,255);

    double norm_factor = 1 / (double) (DIM*DIM);
    #ifdef GIF
        GifWriter gif_writer;
        GifBegin(&gif_writer, gif_filename, DIM, DIM, gif_delay);
        write_gif_frame(lat, &gif_writer, gif_delay);
    #endif
    COLOR colors[2] = {COLOR::white, COLOR::black};
    for (int i=0; i<args.sweeps; i++) {
        for (const auto &color : colors) {
            if (progress != nullptr) {
                progress->write((double)i/args.sweeps);

            }
            broadcast_lattice();

            //cout << "Met: " << action << " vs " << full_action() << endl;
            dS = sweep(color);
            collect_changes(dS, color);
        }
        //cout << "Met: " << action << " vs " << full_action() << endl;

        //if (i==args.thermalization && gif)
            //GifWriteFrame(&gif_writer, white_vec.data(), DIM, DIM, gif_delay); // add one white frame after thermalization

        if (i%args.record_rate==0 && i>=args.thermalization) {
            flow(args.ts, recorder, 0.001);
            #ifdef GIF
                write_gif_frame(lat, &gif_writer, gif_delay);
            #endif
            //if (gif) write_gif_frame(lat, &gif_writer, gif_delay);
        }

        if (i%args.cluster_rate==0 && args.cluster_algorithm == WOLFF) {
            //cout << "Wolff: " << action << " vs " << full_action() << endl;
            #ifdef GIF
                write_gif_frame(lat, &gif_writer, gif_delay);
            #endif
            wolff(); 
            #ifdef GIF
                write_gif_frame(lat, &gif_writer, gif_delay);
            #endif
            //if (gif) write_gif_frame(this, &gif_writer, gif_delay);
            //cout << "Wolff: " << action << " vs " << full_action() << endl;
            
        }

    }
    #ifdef GIF
        GifEnd(&gif_writer);
    #endif

}


double Sweeper::sweep(COLOR color){


    double tot_dS = 0;
    double dS, new_L, old_L, A, r;
    Phi newphi, dphi, phi, backward_nphi_sum, forward_nphi_sum;
    int i; 
    site s;
    vector<site> neighbors;

    array<double, N> zero_array = {};
    Phi zero_phi(zero_array);

    for (int i = 0; i<sites_per_node; i++){
        s = mpi_assignments[color]->at(i);

        phi = lat[s];

        neighbors = lat.neighbor_map[s];
        backward_nphi_sum = lat[neighbors[0]] + lat[neighbors[1]];
        forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

        newphi = new_value(phi);

        old_L = lat.lagrangian( phi, forward_nphi_sum);
        new_L = lat.lagrangian( newphi, forward_nphi_sum);

        dphi = newphi - phi;
        dS = (new_L - old_L) - lat.beta * backward_nphi_sum * dphi;
       
        //cout << old_L << " " << new_L << " " << backward_nphi_sum << endl;
        A = exp(-dS);
        r = randf();

        if (dS < 0 || r <= A) { 
            dphis[i] = dphi;
            tot_dS += dS;
        } else {
            dphis[i] = zero_phi;
        }
    }
    return tot_dS;
}



double randf() {
    return (double)rand() / RAND_MAX;
}

int randint(int n) {
    return rand() % n; // may want to replace this with something better later    
}


double Sweeper::Padd(Phi dphi, Phi phi_b){

    double dS = -lat.beta * (dphi * phi_b); 
    return 1 - exp(-dS); // Schaich Eq. 7.17, promoted for vectors
}

unordered_set <int> Sweeper::generate_cluster(int seed, Phi r, bool accept_all) {
    int s, c, i;
    Phi phi_a, phi_b, dphi;

    double cumsum_dS = 0;

    double Padd_val;
    stack <tuple<int, double>> to_test; // (site, previous r_proj
    unordered_set <int> cluster;

    phi_a = lat[seed];
    
    double r_proj_a = phi_a * r;
    double r_proj_b;
    double proj_sign = sign(r_proj_a);

    //to_test.push(make_tuple(seed, phi_a * numeric_limits<double>::max() )); // site and phi value, using infinity to ensure Padd=1 for first addition
    to_test.push(make_tuple(seed, 0. )); // site and r_projection. r_projection is overridden for seed
    bool first = true;
    while (to_test.size()>0) {
        tie(s, r_proj_a) = to_test.top();
        dphi = -2 * r_proj_a * r;
        to_test.pop();

        if (cluster.find(s)!=cluster.end()) { 
            continue; 
        }

        phi_b = lat[s];
        r_proj_b = phi_b * r;
        if (sign(r_proj_b) == proj_sign) {
            if (accept_all || first || randf() < Padd(dphi, phi_b)) {
                cluster.insert(s);
                for (const int n : lat.neighbor_map[s]){
                    to_test.push( make_tuple(n, r_proj_b) );
                }
                //if (s == seed) {
                    //cout << "SEED" << endl;
                    //cumsum_dS += 2 * beta * (lat[seed] * phi_b);
                //} else { 
                    //cumsum_dS += 2 * beta * (phi_a * phi_b);
                //}
            }

        }
        if (first) first = !first;
    }
    //cout << "cumsum_dS: " << cumsum_dS << endl;
    return cluster;
}

void Sweeper::wolff() {

    int seed = randint(pow(DIM,2));
    unordered_set <int>  cluster;
    int neighbors[4];
    int n, i, c;

    Phi r = random_phi();
    cluster = generate_cluster(seed, r, false);
    
    double dS = 0;
    Phi phi, dphi;


    for (const int c : cluster) {
        phi = lat[c];
        dphi = -2 * (phi * r) * r;
        lat[c] = phi + dphi;
        for (const int n : lat.neighbor_map[c]) {
            dS -= lat.beta * (lat[n] * dphi);
        }
    }

    lat.action+=dS;
    //cout << "dS: " << dS << endl;
#ifdef VERIFY_ACTION
    assert_action();
#endif

}

void Sweeper::broadcast_lattice() {
    if (size_Of_Cluster>1) {
        const int raw_data_len = DIM*DIM*N;
        double raw_data[raw_data_len];

        if (process_Rank == MASTER) {
            for (int i = 0; i < raw_data_len; i++) {
                raw_data[i] = lat[i/N][i%N];
            };
        }


        MPI_Bcast(&raw_data, raw_data_len, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&lat.action, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


        if (process_Rank != MASTER) {
            for (int i = 0; i < raw_data_len; i++) {
                lat[i/N] [i%N] = raw_data[i];
            };
        }
    }

}

void Sweeper::collect_changes(double dS, COLOR color){
    const int recv_data_size = N*DIM*DIM/2;
    double send_data[N*sites_per_node];

    for (int i=0; i<sites_per_node; i++) {
        for (int j = 0; j<N; j++) {
            send_data[i*N+j] = dphis[i][j] ;
        }; 
    };
    
    int recv_sites[DIM*DIM/2];
    double recv_data[recv_data_size];
    double recv_actions[size_Of_Cluster];

    MPI_Gather(mpi_assignments[color]->data(), sites_per_node, MPI_INT, &recv_sites, sites_per_node, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Gather(&send_data, N*sites_per_node, MPI_DOUBLE, &recv_data, N*sites_per_node, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Gather(&dS, 1, MPI_DOUBLE, &recv_actions, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (process_Rank == MASTER) {
        for (int i = 0; i<DIM*DIM/2; i++){
            for (int j = 0; j<N; j++) {
                lat[recv_sites[i]][j] += recv_data[i*N + j];
            }
        }

        for (int i = 0; i<size_Of_Cluster; i++) {
            lat.action += recv_actions[i];
        }
    }
}

Phi sign(Phi x){ // temporary, should be eventually replaced with proj_vec
    Phi phi_sign;
    phi_sign[0] = (x[0] > 0) - (x[0] < 0);
    return phi_sign;
}

double sign(double x){
    return (x>0) - (x<0);
}


// GF
//
//

inline void deriv(Lattice2D& f, double t, const Lattice2D& yn, double h, const Lattice2D* k = nullptr) {

    Phi neighbor_sum; 
    Phi dte;
    Phi e;
    double Pij;
    double laplacianj;
    int L = yn.L;

    for (site s=0; s<L*L; s++) {
        neighbor_sum.init_as_zero();

        if (k) {
            e = yn[s] + k->at(s);
            for (site n : f.neighbor_map[s])
                neighbor_sum += yn[n] + k->at(n);
        } else {
            e = yn[s];
            for (site n : f.neighbor_map[s])
                neighbor_sum += yn[n];
        }

        //if (s==0) cout << "Site: " << e << endl;

        for (int i=0; i<N; i++) {
            dte[i] = 0;
            for (int j=0; j<N; j++) {
                Pij = (i==j) - e[i] * e[j];
                laplacianj = neighbor_sum[j] - 2*D*e[j];
                //if (s==0 && i==0 && j==0) cout << "Laplacian:    " << laplacianj << endl;
                //if (s==0 && i==0 && j==0) cout << "Pij:          " << Pij << endl;
                //if (s==0 && i==0 && j==0) cout << "Neighbor sum: " << neighbor_sum << endl;
                dte[i] += Pij * laplacianj;
            }
        }
        //if (s==0) cout << "deriv " << dte << endl;
        //if (s==0) cout << "deriv " << h*dte << endl;

        f[s] = h*dte;
    }
}

void Sweeper::runge_kutta(double t_, double h, Lattice2D& l, bool recycle_k1) {

        // Runge Kutta (see http://www.foo.be/docs-free/Numerical_Recipe_In_C/c16-1.pdf)
        // Slight changes for efficency:
        //   - deriv(t, y, h, k) := h * f(t, y + k);
        //   - k1 => k1/2; k2 => k2/2
        
        if (t_>0) {
            if (recycle_k1) {
                k1 /= 2;
            } else {
                deriv(k1, t_,     l, h/2);
            }
            deriv(k2, t_+h/2, l, h/2, &k1);
            deriv(k3, t_+h/2, l, h,   &k2);
            deriv(k4, t_+h,   l, h,   &k3);

            k1 /= 3;
            k4 /= 6;

            l += k1;
            l += k2;
            l += k3;
            l += k4;
            
            // Normalize phi
            for (Phi& phi : l) {
                phi /= sqrt(phi.norm_sq());
            }
        }
}

void Sweeper::flow(vector<double> ts, Recorder* recorder, double max_error) {
    // ts must be in ascending order
    double h = 0.01; // aka dt
    double t_ = 0;
    double chi_m;
    double S;

    double error;

    auto measurement_iter = ts.begin();
    double measurement_t = *measurement_iter;

    flowed_lat = lat;

    const auto gif_filename = "flow.gif";
    const int gif_delay = 10;
#ifdef GIF
        GifWriter gif_writer;
        GifBegin(&gif_writer, gif_filename, DIM, DIM, gif_delay);
#endif

#ifdef GIF
        int counter = 0;
#endif

    bool rerun=false;

    while (true) {

        if (t_ + h > measurement_t) {
            runge_kutta(t_, measurement_t-t_, flowed_lat);
            recorder->record(flowed_lat, measurement_t);
            //cout << "took measurement at " << t_ << endl;
            t_ = measurement_t;

            measurement_iter++;
            if (measurement_iter == ts.end()) break; 
            measurement_t = *(measurement_iter);

        } else {
#ifdef ADAPTIVESTEP
            if (!rerun){
                prev_flowed_lat = flowed_lat;
                flowed_lat_2 = flowed_lat;
                runge_kutta(t_, h, flowed_lat);
            } else {
                flowed_lat = flowed_lat_2;
            }

            runge_kutta(t_, h/2, flowed_lat_2, true);
            runge_kutta(t_+h/2, h/2, flowed_lat_2);
            
            error = 0;
            for (site i=0; i<flowed_lat.size(); i++) {
                error += (flowed_lat[i] - flowed_lat_2[i]).norm_sq();
            }

            error = sqrt(error/flowed_lat.size())/15;
            //cout << "t_: "<<t_<< "\terror: " << error << "\th: " << h << endl;
            if (error>max_error) {
                //cout << "RERUN" << endl;
                rerun=true;
                h /= 2;
            } else if (error<max_error/2) {
                rerun=false;
                h *=2;
                t_ += h;
            } else {
                rerun=false;
                t_ += h;
            }

#else
            runge_kutta(t_, h, flowed_lat);
            t_ += h;
#endif
        }
#ifdef ACTIONFLOW
        flowed_lat.action = flowed_lat.full_action();
#endif


#ifdef GIF
            if (counter % 10 == 0) write_gif_frame(flowed_lat, &gif_writer, gif_delay);
            counter++;
#endif
    }

    #ifdef GIF
        GifEnd(&gif_writer);
        system("gifsicle --colors 256 --resize 512x512 flow.gif -o flow.gif");
    #endif
}

