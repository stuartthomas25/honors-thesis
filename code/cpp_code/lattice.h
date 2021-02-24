#ifndef LATTICE_H
#define LATTICE_H


#include "phi.h"
#include <vector>

typedef tuple<int,int,int,int> Plaquette;

using namespace std;

class Lattice2D {
    typedef vector<Phi> datatype;
    datatype data;
    public:
        static int L;
        static double beta;
        static vector<vector<int>> neighbor_map;
        static vector<Plaquette> plaquette_map;
        double action;

        vector<Phi> vec() const;
        size_t size() const;

        Lattice2D();
        Lattice2D(const Lattice2D& other);

        Phi operator[] (int i) const;
        Phi& operator[] (int i);
        Phi at(int i) const;

        Lattice2D& operator+=(const Lattice2D& other);
        Lattice2D& operator*=(const double & factor);
        Lattice2D& operator/=(const double & factor);
        Lattice2D operator+ (const Lattice2D & other) const;
        //Lattice2D operator+ (const Lattice2D & other) const;

        // Iterator stuff
        typedef datatype::iterator iterator;
        typedef datatype::const_iterator const_iterator;

        iterator begin();
        const_iterator cbegin() const;
        iterator end();
        const_iterator cend() const;

        static double lagrangian(const Phi phi, const Phi nphi_sum);
        double full_action();

};
#endif
