#include <iostream>
#include <stdexcept>
#include "sweep.h"
#include <algorithm>
#include "gif.h"
#include "progress.cpp"
#include <assert.h>
#include <stdlib.h>

#define VERIFY_ACTION 0

using namespace std;

typedef int site;

double Sweeper::beta;





































































































Sweeper::Sweeper (int DIM, MPI_Comm c, bool makeGif) : 
    process_Rank(get_rank(c)), 
    size_Of_Cluster(get_size(c)), 
    sites_per_node((DIM*DIM)/2 / get_size(c)),
    DIM{DIM},
    gif{makeGif}
{
    // generate lattice
    array<double, N> zero_array = {};
    Phi zero_phi(zero_array);

    int site;
    for (int i = 0; i<DIM; i++){
        for (int j = 0; j<DIM; j++) {
            site = i*DIM + j;
            lat[site] = new_value(zero_phi);
            
            lat.neighbor_map[site] = full_neighbors(site); 

            // generate MPI assignments
            if (site / (2*sites_per_node) == process_Rank) { // Factor of 2 for checkerboard
                if ( (i+j)%2 ) {
                    mpi_assignments[COLOR::black].push_back(site);
                } else {
                    mpi_assignments[COLOR::white].push_back(site);
                }   
            }
        }
    }

    action = full_action(); 

    int offset = process_Rank *  sites_per_node;

}

Sweeper::~Sweeper() = default;

void print(std::vector <Phi> const &a) {
   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
   cout << endl;
}

double Sweeper::full_action() {
    int site, i;
    Phi forward_nphi_sum;
    vector<int> neighbors;
    double S = 0;
    // initial action
    for (int s=0; s<DIM*DIM; s++) {
        neighbors = lat.neighbor_map[s];
        forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

        S += lagrangian(lat[s], forward_nphi_sum);
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

vector<int> Sweeper::full_neighbors(int site) {

    int x = site % DIM;
    int y = site / DIM;
    vector<int> neighbors; 
    neighbors.push_back(wrap(x-1) + DIM * y);
    neighbors.push_back(x + DIM * wrap(y-1));
    neighbors.push_back(wrap(x+1) + DIM * y);
    neighbors.push_back(x + DIM * wrap(y+1));
    return neighbors;
}


double Sweeper::lagrangian(Phi phi, Phi nphi_sum) {
    return beta * (D - phi * nphi_sum); // note that the sum over dimension has already been made
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
    double fa;
    if (VERIFY_ACTION) {
        fa = full_action();
        if (abs(fa-action)>tol) {
            cout << "ASSERT ACTION FAILED (" << action << " != " << fa << ")\n";
            exit(1);
        } else {
            cout << "assert action passed (" << action << " == " << fa << ")\n";

        }
    }
}

vector<measurement> Sweeper::full_sweep(const sweep_args& args = sweep_args()) {
   
    if (args.cluster_algorithm != WOLFF) throw invalid_argument("Currently, Wolff is the only allowed algorithm");

    shared_ptr<progress_bar> progress;
    progress = args.progress ? make_shared<progress_bar>(cout, 70u, "Working") :  nullptr; 

    vector<measurement> measurements;
    Lattice2D dphis;
    double dS;
    Phi phibar;
    int s;

    GifWriter gif_writer;
    if (gif) GifBegin(&gif_writer, gif_filename, DIM, DIM, gif_delay);

    vector<uint8_t> white_vec(DIM*DIM*4,255);

    double norm_factor = 1 / (double) (DIM*DIM);
    if (gif) write_gif_frame(lat, &gif_writer, gif_delay);
    COLOR colors[2] = {COLOR::white, COLOR::black};
    for (int i=0; i<args.sweeps; i++) {
        for (const auto &color : colors) {
            if (progress != nullptr) {
                progress->write((double)i/args.sweeps);

            }
            broadcast_lattice();

            //cout << "Met: " << action << " vs " << full_action() << endl;
            auto [dphis, dS] = sweep(color);
            collect_changes(dphis, dS, color);
            //assert_action();
        }
        //cout << "Met: " << action << " vs " << full_action() << endl;

        //if (i==args.thermalization && gif)
            //GifWriteFrame(&gif_writer, white_vec.data(), DIM, DIM, gif_delay); // add one white frame after thermalization

        if (i%args.record_rate==0 && i>=args.thermalization) {
            phibar.init_as_zero();
            for (const Phi& phi : lat) {
                phibar = phibar + phi;
            }
            measurement aMeasurement = {phibar*norm_factor, action};
            measurements.push_back(aMeasurement);
            if (gif) write_gif_frame(lat, &gif_writer, gif_delay);
            
        }

        if (i%args.cluster_rate==0) {
            //cout << "Wolff: " << action << " vs " << full_action() << endl;
            wolff(); 
            //assert_action();
            //if (gif) write_gif_frame(this, &gif_writer, gif_delay);
            //cout << "Wolff: " << action << " vs " << full_action() << endl;
            
        }

    }
    if (gif) GifEnd(&gif_writer);
    return measurements;

}


tuple<vector<Phi>, double> Sweeper::sweep(COLOR color){

    vector<Phi> dphis(sites_per_node);

    double tot_dS = 0;
    double dS, new_L, old_L, A, r;
    Phi newphi, dphi, phi, backward_nphi_sum, forward_nphi_sum;
    int i, site;
    vector<int> neighbors;

    array<double, N> zero_array = {};
    Phi zero_phi(zero_array);

    for (int i = 0; i<sites_per_node; i++){
        site = mpi_assignments[color] [i];

        phi = lat[site];

        neighbors = lat.neighbor_map[site];
        backward_nphi_sum = lat[neighbors[0]] + lat[neighbors[1]];
        forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

        newphi = new_value(phi);
        if (abs(newphi[0])>1e+10){
            exit(1);
        }

        old_L = lagrangian( phi, forward_nphi_sum);
        new_L = lagrangian( newphi, forward_nphi_sum);

        dphi = newphi - phi;
        dS = (new_L - old_L) - beta * backward_nphi_sum * dphi;
       
        //cout << old_L << " " << new_L << " " << backward_nphi_sum << endl;
        A = exp(-dS);
        r = randf();

        if (dS < 0 || r <= A) { 
        //if (dS < 0 ) { 
            dphis[i] = dphi;
            tot_dS += dS;
        } else {
            dphis[i] = zero_phi;
        }
    }
    return make_tuple(dphis, tot_dS);

}



double randf() {
    return (double)rand() / RAND_MAX;
}

int randint(int n) {
    return rand() % n; // may want to replace this with something better later    
}

double Sweeper::Padd(Phi phi_a, Phi phi_b){
    return 1 - exp(-2*(phi_a*phi_b)); // Schaich Eq. 7.17, promoted for vectors
}

unordered_set <int> Sweeper::generate_cluster(int seed, bool accept_all) {
    int s, c, i;
    Phi phi_a, phi_b;

    double Padd_val;
    stack <tuple<int, Phi>> to_test;
    unordered_set <int> cluster;

    phi_a = lat[seed];
    
    Phi r = random_phi();
    double r_projection = phi_a * r;
    double proj_sign = sign(r_projection);
    to_test.push(make_tuple(seed, phi_a * numeric_limits<double>::max() )); // site and phi value, using infinity to ensure Padd=1 for first addition
    while (to_test.size()>0) {
        tie(s, phi_a) = to_test.top();
        to_test.pop();

        if (cluster.find(s)!=cluster.end()) { 
            continue; 
        }

        phi_b = lat[s];
        if (sign(r * phi_b) == proj_sign) {
            Padd_val = Padd(phi_a, phi_b);
            if (accept_all || randf() < Padd_val) {
                cluster.insert(s);
                for (const int n : lat.neighbor_map[s]){
                    to_test.push( make_tuple(n, phi_b) );
                }
            }

        }
    }
    return cluster;
}

void Sweeper::wolff() {

    int seed = randint(pow(DIM,2));
    unordered_set <int>  cluster;
    int neighbors[4];
    int n, i, c;

    cluster = this->generate_cluster(seed, false);
    
    double dS = 0;
    Phi phi;


    unordered_set<int>::iterator it;
    for (it = cluster.begin(); it != cluster.end(); ++it) {
        c = *it;
        phi = lat[c];
        lat[c] = -phi;
        for (const int n : lat.neighbor_map[c]) {
            dS += 2 * beta * (lat[n] * phi);
        }
    }

    action+=dS;

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
        MPI_Bcast(&action, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


        if (process_Rank != MASTER) {
            for (int i = 0; i < raw_data_len; i++) {
                lat[i/N] [i%N] = raw_data[i];
            };
        }
    }

}

void Sweeper::collect_changes(vector<Phi> dphis, double dS, COLOR color){
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

    MPI_Gather(mpi_assignments[color].data(), sites_per_node, MPI_INT, &recv_sites, sites_per_node, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Gather(&send_data, N*sites_per_node, MPI_DOUBLE, &recv_data, N*sites_per_node, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Gather(&dS, 1, MPI_DOUBLE, &recv_actions, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (process_Rank == MASTER) {
        for (int i = 0; i<DIM*DIM/2; i++){
            for (int j = 0; j<N; j++) {
                lat[recv_sites[i]][j] += recv_data[i*N + j];
            }
        }

        for (int i = 0; i<size_Of_Cluster; i++) {
            action += recv_actions[i];
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

void deriv(Lattice2D& f, double t, const Lattice2D& yn, double h, const Lattice2D* k = nullptr) {

    Phi neighbor_sum; 
    Phi dte;
    Phi e;
    double Pij;
    double laplacianj;
    int L = yn.L;

    for (site s=0; s<L*L; s++) {
        neighbor_sum.init_as_zero();

        if (k) {
            e = yn[s] + (*k)[s];
            for (site n : f.neighbor_map[s])
                neighbor_sum += yn[n] + (*k)[n];
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

void Sweeper::flow(double t) {
    double h = 0.01; // aka dt
    double h_2 = h/2;
    double t_ = 0;
    Lattice2D flowed_lat(lat);
    Lattice2D k1, k2, k3, k4;

    const auto gif_filename = "flow.gif";
    const int gif_delay = 10;
    GifWriter gif_writer;
    GifBegin(&gif_writer, gif_filename, DIM, DIM, gif_delay);

    int counter = 0;

    while (t_<t) {
        // Runge Kutta (see http://www.foo.be/docs-free/Numerical_Recipe_In_C/c16-1.pdf)
        // Slight changes for efficency:
        //   - deriv(t, y, h, k) := h * f(t, y + k);
        //   - k1 => k1/2; k2 => k2/2
        
        deriv(k1, t_,     flowed_lat, h_2);
        deriv(k2, t_+h_2, flowed_lat, h_2, &k1);
        deriv(k3, t_+h_2, flowed_lat, h,   &k2);
        deriv(k4, t_+h,   flowed_lat, h,   &k3);

        k1 /= 3;
        k2 /= (3/2);
        k3 /= 3;
        k4 /= 6;

        flowed_lat += k1;
        flowed_lat += k2;
        flowed_lat += k3;
        flowed_lat += k4;

        //cout << flowed_lat[0] << endl;
        //for (auto n : flowed_lat.neighbor_map[0]) {
            //cout << flowed_lat[n] << "\t";
        //}
        //cout << endl;
        cout << t_ << " " << flowed_lat.full_action(Sweeper::lagrangian) << endl;
        //cout << endl;

        //deriv(k1, t_,     flowed_lat, h);

        //cout << flowed_lat[0] << endl;
        //cout << k1[0] << endl;
        //flowed_lat += k1;
        //cout << flowed_lat[0] << endl;
        //cout << endl;
        
        // Normalize phi
        for (Phi& phi : flowed_lat) {
            phi /= sqrt(phi.norm_sq());
        }

        t_ += h;

        if (counter % 10 == 0) write_gif_frame(flowed_lat, &gif_writer, gif_delay);
        counter++;
    }

    GifEnd(&gif_writer);
    system("gifsicle --colors 256 --resize 512x512 flow.gif -o bigger.gif");

}

