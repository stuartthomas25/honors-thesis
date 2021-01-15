#include <iostream>
#include "sweep.h"
#include <algorithm>

using namespace std;

namespace croutines {

    Phi::Phi() {}; 

    Phi::Phi(array<double, N> phi){
        this->phi = phi;
    };

    void Phi::init_as_zero(){
        for (int i=0; i<N; i++) {
            phi[i] = 0;
        }
    }

    double Phi::norm_sq() const {
        double cumsum = 0;
        for (int i=0; i<N; i++) {
            cumsum += phi[i] * phi[i];
        }
        return cumsum;

    }


    Phi Phi::operator+ (const Phi & aPhi) const {
        Phi new_phi;
        for (int i=0; i<N; i++) {
            new_phi[i] = phi[i] + aPhi[i];
        }
        return new_phi;
    }

    Phi Phi::operator- (const Phi & aPhi) const {
        return *this + (-aPhi);    
    }
    Phi Phi::operator- () const {
        Phi new_phi;
        for (int i=0; i<N; i++) {
            new_phi[i] = -phi[i];
        }
        return new_phi;
    }
    double Phi::operator* (const Phi & aPhi) const{
        //dot product
        double dot = 0;
        for (int i=0; i<N; i++) {
            dot += phi[i] * aPhi[i];
        }
        return dot;
    }

    Phi Phi::operator* (const double & a) const{
        Phi new_phi;
        for (int i=0; i<N; i++) {
            new_phi[i] = phi[i] * a;
        }
        return new_phi;

    }
    void Phi::print(){
    }

    ostream& operator<<(ostream& os, const Phi & aPhi) 
    { 
        os << "(";
        for (int i=0; i<N; i++) {
            os << aPhi[i];
            if (i<N-1) { os<<", "; }
        }
        os << ")";
        return os; 
    }
    
    double Phi::operator[] (int i) const {
        return phi[i];
    }
    double & Phi::operator[] (int i) {
        return phi[i];
    }

    bool Phi::operator== (const Phi & aPhi) const {
        for (int i=0; i<N; i++) {
            if (phi[i] != aPhi[i]) {
                return false;
            }
        }
        return true;
    }


    Sweeper::Sweeper (double m02, double lam, int DIM, MPI_Comm c, bool gif=false) : 
        process_Rank(get_rank(c)), 
        size_Of_Cluster(get_size(c)), 
        sites_per_node((DIM*DIM)/2 / get_size(c)),
        DIM{DIM},
        gif{gif}

    {

        
        this->redef_mass = 2 + 0.5 * m02;
        this->quarter_lam = 0.25 * lam;
        this->beta=1; //EDIT


        if (gif) {
            GifBegin(&gif_writer, gif_filename, DIM, DIM, gif_delay);
        }

        // generate lattice
        array<double, N> zero_array = {};
        Phi zero_phi(zero_array);

        for (int i=0; i<DIM*DIM; i++) {
            lat.push_back(new_value(zero_phi));
            if (gif) {
                for (int j=0; j<4; j++) {
                    gif_vec.push_back(0);
                }
            }
        }


        int site;

        action = full_action(); 
        // generate MPI assignments
        for (int i = 0; i<DIM; i++){
            for (int j = 0; j<DIM; j++) {
                site = i*DIM + j;
                if ( (i+j)%2 ) {
                    mpi_assignments[COLOR::black].push_back(site);
                } else {
                    mpi_assignments[COLOR::white].push_back(site);
                }   
            }
        }
       
        int offset = process_Rank *  sites_per_node;

    }

    Sweeper::~Sweeper(){
        GifEnd(&gif_writer);
    }

    double Sweeper::full_action() {
        int site, i;
        Phi forward_nphi_sum;
        int neighbors[4];
        double S = 0;
        // initial action
        for (int s=0; s<DIM*DIM; s++) {
            full_neighbors(s, neighbors);
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

    void Sweeper::full_neighbors(int site, int neighbors[4]) {

        int x = site % DIM;
        int y = site / DIM;
        neighbors[0] = wrap(x+1) + DIM * y;
        neighbors[1] = wrap(x-1) + DIM * y;
        neighbors[2] = x + DIM * wrap(y+1);
        neighbors[3] = x + DIM * wrap(y-1);
    }


    double Sweeper::lagrangian(Phi phi, Phi nphi_sum) {
        return this->phi4_lagrangian(phi, nphi_sum);    
    }
    double Sweeper::phi4_lagrangian(Phi phi, Phi nphi_sum) {
        double phi2 = phi*phi;
        return -phi * nphi_sum + this->redef_mass * phi2 + this->quarter_lam * phi2*phi2;
    }

    double Sweeper::sigma_lagrangian(Phi phi, Phi nphi_sum) {
        return -1 * this->beta * (phi * nphi_sum);
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
        if (N!=1) { throw invalid_argument("This should only be used for the Phi4 model"); }

        Phi dphi;
        dphi[0] = rand_dist(randf());
        return dphi;
    }

    Phi Sweeper::proj_vec () {
        //for (int i = 0; i<N-1; i++) {
            //r[i] = 
        //}
        return Phi();
    }

    Phi Sweeper::random_phi() {
        bool phi4 = true;
        if (phi4) {
            Phi ret_phi({3 * randf() - 1.5});
            return ret_phi;
        } else {
            Phi new_phi;
            for (int i=0; i<N; i++) {
                new_phi[i] = 2 * randf() - 1;
            }
            return new_phi * (1/sqrt(new_phi.norm_sq()));
        }  

    }

    void Sweeper::write_gif_frame(double rate=3) {
        int val;
        Phi aPhi;
        array<uint8_t, 4> pixel = {0,0,0,0};
        fill(gif_vec.begin(), gif_vec.end(), 50);

        if (gif) {
            for (int i=0; i<DIM*DIM; i++) {
                aPhi = lat[i];
                val = 255 * (1 - exp(-rate * aPhi[0])) / (1 + exp(-rate * aPhi[0]));
                if (val >= 0) {
                    gif_vec[4*i] = static_cast<uint8_t>(val);
                } else {
                    gif_vec[4*i+1] = static_cast<uint8_t>(-val);
                }

                gif_vec[4*i + 3] = 0; 
                
            }
            GifWriteFrame(&gif_writer, gif_vec.data(), DIM, DIM, gif_delay);
        }
    }


    vector<measurement> Sweeper::full_sweep(const sweep_args& args = sweep_args()) {
       
        if (args.cluster_algorithm != WOLFF) { throw invalid_argument("Currently, Wolff is the only allowed algorithm"); }

        shared_ptr<progress_bar> progress;
        if (args.progress) {
            progress = make_shared<progress_bar>(cout, 70u, "Working");
        } else {
            progress = nullptr; 
        }

        vector<measurement> measurements;

        vector<Phi> dphis;
        double dS;
        Phi phibar;
        int s;

        vector<uint8_t> white_vec(DIM*DIM*4,255);

        double norm_factor = 1 / (double) (DIM*DIM);
        write_gif_frame();
        COLOR colors[2] = {COLOR::white, COLOR::black};
        for (int i=0; i<args.sweeps; i++) {
            for (const auto &color : colors) {

                if (progress != nullptr) {
                    progress->write((double)i/args.sweeps);

                }
                broadcast_lattice();

                //cout << "Met: " << action << " vs " << full_action() << endl;
                tie(dphis, dS) = sweep(color);
                collect_changes(dphis, dS, color);
            }
            write_gif_frame();
            //cout << "Met: " << action << " vs " << full_action() << endl;

            if (i==args.thermalization && gif){
                GifWriteFrame(&gif_writer, white_vec.data(), DIM, DIM, gif_delay);
            }

            if (i%args.record_rate==0 && i>=args.thermalization) {
                phibar.init_as_zero();
                for (Phi l : lat) {
                    phibar = phibar + l;
                }
                measurement aMeasurement = {phibar*norm_factor, action};
                measurements.push_back(aMeasurement);
                
            }

            if (i%args.cluster_rate==0) {
                //cout << "Wolff: " << action << " vs " << full_action() << endl;
                wolff();
                write_gif_frame();
                //cout << "Wolff: " << action << " vs " << full_action() << endl;
                
            }

        }
        return measurements;

    }


    tuple<vector<Phi>, double> Sweeper::sweep(COLOR color){

        vector<Phi> dphis;

        double tot_dS = 0;
        double dS, new_L, old_L, A, r;
        Phi dphi, phi, backward_nphi_sum, forward_nphi_sum;
        int i, site;
        int neighbors[4];

        array<double, N> zero_array = {};
        Phi zero_phi(zero_array);

        for (int i = 0; i<sites_per_node; i++){
            site = mpi_assignments[color] [i];

            phi = this->lat[site];

            this->full_neighbors(site, neighbors);
            backward_nphi_sum = lat[neighbors[0]] + lat[neighbors[1]];
            forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

            dphi = this->new_value(phi);
            old_L = this->lagrangian( phi, forward_nphi_sum);
            new_L = this->lagrangian( phi+dphi, forward_nphi_sum);
            dS = (new_L - old_L) - backward_nphi_sum * dphi;
           
            //cout << old_L << " " << new_L << " " << backward_nphi_sum << endl;
            A = exp(-dS);
            r = randf();

            if (dS < 0 || r <= A) { 
            //if (dS < 0 ) { 
                dphis.push_back(dphi);
                tot_dS += dS;
            } else {
                dphis.push_back(zero_phi);
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

    set <int> Sweeper::generate_cluster(int seed, bool accept_all) {
        int s, c, i;
        Phi phi_a, phi_b;

        double Padd_val;
        stack <tuple<int, Phi>> to_test;
        int neighbors[4];
        set <int> cluster;

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
                    full_neighbors(s, neighbors);
                    for (i=0; i<4; i++){
                        to_test.push( make_tuple(neighbors[i], phi_b) );
                    }
                }

            }
        }
        return cluster;
    }

    //tuple<vector<Phi>, double> Sweeper::wolff() {
    void Sweeper::wolff() {

        int seed = randint(pow(DIM,2));
        set <int>  cluster;
        int neighbors[4];
        int n, i, c;

        cluster = this->generate_cluster(seed, false);
        
        double dS = 0;
        Phi phi;


        set<int>::iterator it;
        for (it = cluster.begin(); it != cluster.end(); ++it) {
            c = *it;
            phi = lat[c];
            lat[c] = -phi;
            this->full_neighbors(c, neighbors);
            for (i=0; i<4; i++) {
                n = neighbors[i];
                dS += 2 * (lat[n] * phi);
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

        for (int i = 0; i<sites_per_node; i++) {
            for (int j = 0; j<N; j++) {
                send_data[i*N+j] = dphis[i][j] ;
            }; 
        };
        
        int recv_sites[DIM*DIM/2];
        double recv_data[recv_data_size];
        double recv_actions[size_Of_Cluster];



        //cout << 0;
        MPI_Gather(mpi_assignments[color].data(), sites_per_node, MPI_INT, &recv_sites, DIM*DIM/2, MPI_INT, MASTER, MPI_COMM_WORLD);
        //cout << 1;
        MPI_Gather(&send_data, N*sites_per_node, MPI_DOUBLE, &recv_data, recv_data_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        //cout << 2;
        MPI_Gather(&dS, 1, MPI_DOUBLE, &recv_actions, size_Of_Cluster, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        //cout << 3;

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

    /*double dot(array<double, N> a, array<double, N> b){*/
        //double sum = 0;
        //for (int i = 0; i<N; i++) {
            //sum += a[i] * b[i];
        //}
    /*}*/


}
