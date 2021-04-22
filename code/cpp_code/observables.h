#include <math.h>
#include <string>


namespace observables {


    class action : public BaseObservable {
        public:
            double operator()(const Lattice2D& lat) const {
                return lat.action;
            };
            const string name() const { return "S"; };
    };

    class beta : public BaseObservable {
        public:
            double operator()(const Lattice2D& lat) const {
                return lat.beta;
            };
            const string name() const { return "beta"; };
    };

    class L : public BaseObservable {
        public:
            double operator()(const Lattice2D& lat) const {
                return lat.L;
            };
            const string name() const { return "L"; };
    };

    class chi_m : public BaseObservable {
        public:
            double operator()(const Lattice2D& lat) const {
                double val = 0;

                for (auto itx = lat.cbegin(); itx!=lat.cend(); ++itx) {
                    for (auto ity = lat.cbegin(); ity!=lat.cend(); ++ity) {
                        val += (*itx)*(*ity);
                    }
                }
                return val;
            };
            const string name() const { return "chi_m"; };
    };

    class F : public BaseObservable {
        public:
            double operator()(const Lattice2D& lat) const {
                double val = 0;

                int x=0;
                int y=0;
                const double two_pi_L = 2*M_PI*lat.L;

                //for (const Phi& phi : as_const(lat)) {
                for (auto itx = lat.cbegin(); itx!=lat.cend(); ++itx) {
                    for (auto ity = lat.cbegin(); ity!=lat.cend(); ++ity) {
                        val += (*itx)*(*ity) * cos( two_pi_L * (x - y));
                        x++;
                        if (x==lat.L) x=0; // This ensures that x is the Euclidean space dimension
                    }
                    if (y==lat.L) y=0;
                }
                return val;
            };
            const string name() const { return "F"; };
    };

    class Q : public BaseObservable {
        private:
             double angle(double re, double im) const {
                double arctan = atan(im / re);
                if (re> 0) {
                    return arctan;
                } else if (im> 0) {
                    return M_PI + arctan;
                } else {
                    return (-M_PI + arctan);
                }
            }
            double sigma_A(Phi s1, Phi s2, Phi s3) const {
                // Returns values (-2pi,2pi)
                double real_part = 1 + s1 * s2 + s2 * s3 + s3 * s1;
                double imag_part = s1 * (s2 & s3);
                return 2*angle(real_part, imag_part);
            }

            double q(int x, const Lattice2D& lat, bool reversed=false) const {
                int x1, x2, x3, x4;
                tie(x1, x2, x3, x4) = lat.plaquette_map[x]; 
                Phi s1 = lat[x1], s2 = lat[x2], s3 = lat[x3], s4 = lat[x4];
                if (reversed) {
                    return sigma_A(s1,s2,s4) + sigma_A(s2,s3,s4);
                } else {
                    return sigma_A(s1,s2,s3) + sigma_A(s1,s3,s4);
                }
            }
        public:
            double operator()(const Lattice2D& lat) const{
                double Q = 0;
                for (int i=0; i<lat.L*lat.L; i++) {
                    Q += q(i, lat);
                }
                Q /= (4 * M_PI);
                return Q;
            };
            const string name() const { return "Q"; };
    };

};
