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
    };



    Phi Phi::operator+ (const Phi & aPhi) const {
        Phi new_phi;
        for (int i=0; i<N; i++) {
            new_phi[i] = phi[i] + aPhi[i];
        }
    }

    Phi Phi::operator- (const Phi & aPhi) const {
        return *this + (-aPhi);    
    }
    Phi Phi::operator- () const {
        Phi new_phi;
        for (int i=0; i<N; i++) {
            new_phi[i] = -phi[i];
        }
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
        cout << "[";
        for (int i=0; i<N; i++) {
            cout << phi[i];
            if (i<N-1) { cout<<", "; }
        }
        cout << "]\n";
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



    Sweeper::Sweeper () {}

    Sweeper::Sweeper (vector<Phi> lattice, double m02, double lam) {
        this->lat = lattice;
        this->dim = (int) sqrt(lattice.size());
        this->redef_mass = 2 + 0.5 * m02;
        this->quarter_lam = 0.25 * lam;
        this->beta=1; //EDIT
    }

    int Sweeper::wrap( int c) {
        int mod = c % dim;
        if (mod < 0) {
            mod += dim;
        }
        return mod;

    }

    void Sweeper::full_neighbors(int site, int neighbors[4]) {

        int dim = this->dim;
        int x = site % dim;
        int y = site / dim;
        neighbors[0] = wrap(x+1) + dim * y;
        neighbors[1] = wrap(x-1) + dim * y;
        neighbors[2] = x + dim * wrap(y+1);
        neighbors[3] = x + dim * wrap(y-1);
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

    vector<measurement> Sweeper::full_sweep(int sweeps,
                                    int thermalization,
                                    int record_rate,
                                    ClusterAlgorithm cluster_algorithm = WOLFF, 
                                    int cluster_rate = 5
                                ) {
       

        if (cluster_algorithm != WOLFF) { throw invalid_argument("Currently, Wolff is the only allowed algorithm"); }


        progress_bar progress{cout, 70u, "Working"};

        vector<int> sites;
        vector<measurement> measurements;

        for (int i=0; i<lat.size(); i++) {
            sites.push_back(i);
        }
        
        vector<Phi*> dphis;
        double dS;
        double action = 0;
        Phi phibar;

        for (int i=0; i<sweeps; i++) {
            progress.write(i/sweeps);
            tie(dphis, dS) = this->sweep(sites);

            for (int s : sites) {
                lat[s] = lat[s] + (*dphis[s]);
            }
            action += dS;
            if (i%cluster_rate==0) {
                tie(lat, dS) = this->wolff();
                action += dS;
            }
            if (i%record_rate==0 && i>=thermalization) {
                phibar.init_as_zero();
                for (int s : sites) {
                    phibar = phibar + lat[s];
                }
                measurement aMeasurement = {phibar, action};
                measurements.push_back(aMeasurement);
                
            }

        }
        cout << "\n\n";

        return measurements;

    }






    tuple<vector<Phi*>, double> Sweeper::sweep(vector<int>& sites){

        vector<Phi*> dphis;

        double tot_dS = 0;
        double dS, new_L, old_L, A, r;
        Phi dphi, phi, backward_nphi_sum, forward_nphi_sum;
        int i;
        int neighbors[4];

        for (int site : sites){
                
            phi = this->lat[site];

            this->full_neighbors(site, neighbors);
            backward_nphi_sum = lat[neighbors[0]] + lat[neighbors[1]];
            forward_nphi_sum  = lat[neighbors[2]] + lat[neighbors[3]];

            dphi = this->new_value(phi);
            old_L = this->lagrangian( phi, forward_nphi_sum);
            new_L = this->lagrangian( phi+dphi, forward_nphi_sum);
            dS = (new_L - old_L) - backward_nphi_sum * dphi;
            
            A = exp(-dS);
            r = randf();

            if (dS < 0 || r <= A) {
                dphis.push_back(&dphi);
                tot_dS += dS;
            } else {
                dphis.push_back(NULL);
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
        return 1 - exp(-2*(phi_a*phi_b)); // Schaich Eq. 7.17
    }

    set <int> Sweeper::generate_cluster(int seed, bool accept_all) {
        int s, n, c, i;
        Phi phi_sign, phi_a, phi_b;

        double Padd;
        stack <tuple<int, Phi>> to_test;
        int neighbors[4];
        set <int> cluster;

        phi_a = lat[seed];
        
        Phi r;
        phi_sign = sign(phi_a); // get the sign of phi_a
        to_test.push(make_tuple(seed, phi_sign * numeric_limits<double>::max() )); // site and phi value, using infinity to ensure Padd=1 for first addition
        while (to_test.size()>0) {
            tie(s, phi_a) = to_test.top();
            to_test.pop();
            
            if (cluster.find(s)==cluster.end()) { continue; }

            phi_b = lat[s];
            if (sign(phi_b) == phi_sign) {
                Padd = this->Padd(phi_a, phi_b);
                if (accept_all || randf() < Padd) {
                    cluster.insert(s);
                    this->full_neighbors(s, neighbors);
                    for (i=0; i<4; i++){
                        to_test.push( make_tuple(n, phi_b) );
                    }
                }
            }
        }
        return cluster;
    }

    tuple<vector<Phi>, double> Sweeper::wolff() {

        int seed = randint(pow(this->dim,2));
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

        return make_tuple(lat, dS);
    }


    Phi sign(Phi x){ // temporary, should be eventually replaced with proj_vec
        Phi phi_sign;
        phi_sign[0] = (x[0] > 0) - (x[0] < 0);
        return phi_sign;
    }

    /*double dot(array<double, N> a, array<double, N> b){*/
        //double sum = 0;
        //for (int i = 0; i<N; i++) {
            //sum += a[i] * b[i];
        //}
    /*}*/


}
