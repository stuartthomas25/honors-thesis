#ifndef PHI_H // include guard
#define PHI_H

#include "constants.h"
#include <array>
#include <ostream>

using namespace std;

class Phi {
    array<double, N> phi;

    public:
        Phi();
        Phi(array<double, N> phi);
        void init_as_zero();
        double norm_sq() const;
        Phi& operator+=(const Phi& other);
        Phi& operator*=(const double & a);
        Phi& operator/=(const double & a);

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

Phi operator*(double a, const Phi& b);
  

#endif
