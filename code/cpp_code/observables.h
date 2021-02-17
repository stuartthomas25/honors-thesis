namespace observables {
    double action(Lattice2D& lat) {
        return lat.action;
    }


    double chi_m(Lattice2D& lat) {
        double val = 0;
        for (const Phi& phi : lat) {
            for (const Phi& phi2 : lat) {
                val += phi2*phi;
            }
        }
        return val;
    }

    namespace {
        inline double angle(double re, double im) {
            double arctan = atan(im / re);
            if (re> 0) {
                //cout << "arctan 1: " << arctan << endl;
                return arctan;
            } else if (im> 0) {
                //cout << "arctan 2: " << arctan << endl;
                return M_PI + arctan;
            } else {
                //cout << "arctan 3: " << arctan << endl;
                return (-M_PI + arctan);
            }
        }
        inline double sigma_A(Phi s1, Phi s2, Phi s3) {
            // Returns values (-2pi,2pi)
            double real_part = 1 + s1 * s2 + s2 * s3 + s3 * s1;
            double imag_part = s1 * (s2 & s3);
            //cout << 2*angle(real_part, imag_part) << endl;
            assert(abs(2*angle(real_part, imag_part)) < 2*M_PI);
            return 2*angle(real_part, imag_part);
        }

        inline double q(int x, Lattice2D& lat) {
            int x1, x2, x3, x4;
            tie(x1, x2, x3, x4) = lat.plaquette_map[x]; 
            Phi s1 = lat[x1], s2 = lat[x2], s3 = lat[x3], s4 = lat[x4];
            //cout << "q: " <<0.25 * M_PI * (sigma_A(s1,s2,s4) + sigma_A(s2,s3,s4)) << endl;
            //cout << sigma_A(s1,s2,s4) << "  " << sigma_A(s2,s3,s4) << endl;
            return 0.25 * M_PI * (sigma_A(s1,s2,s3) + sigma_A(s1,s3,s4));
            //return 0.25 * M_PI * (sigma_A(s1,s2,s4) + sigma_A(s2,s3,s4));
        }
    }

    double Q(Lattice2D& lat) {
        double Q = 0;
        for (int i=0; i<lat.L*lat.L; i++) {
            Q += q(i, lat);
        }
        //cout << "QQ: " << Q << "\n----------\n\n";
        return Q;
    }
}
