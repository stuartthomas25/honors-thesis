// lattice.cpp

#include <iostream>
#include "lattice.h"

int Lattice2D::L;
unordered_map<int,vector<int>> Lattice2D::neighbor_map;
unordered_map<int, Plaquette> Lattice2D::plaquette_map;

Lattice2D::Lattice2D() {
    data.resize(L*L);
}
    
Lattice2D::Lattice2D(const Lattice2D& other) {
    data = other.vec();
}

vector<Phi> Lattice2D::vec() const { return data; }

Phi Lattice2D::operator[] (int i) const { return data[i]; };
Phi& Lattice2D::operator[] (int i) { return data[i]; };

Lattice2D& Lattice2D::operator+=(const Lattice2D& other) {
    auto iter = other.cbegin(); 
    for_each(data.begin(), data.end(), [&iter](Phi &phi){phi += *(iter++); } );
    return *this;
};

Lattice2D& Lattice2D::operator*=(const double & factor) {
    for_each(data.begin(), data.end(), [factor](Phi &phi){phi *= factor; } );
    return *this;
};

Lattice2D& Lattice2D::operator/=(const double & factor) {
    for_each(data.begin(), data.end(), [factor](Phi &phi){phi /= factor; } );
    return *this;
};

Lattice2D Lattice2D::operator+ (const Lattice2D & other) const {
    Lattice2D new_lat(*this);
    new_lat += other;
    return new_lat;
};
        

Lattice2D::iterator Lattice2D::begin() { return data.begin(); }
Lattice2D::const_iterator Lattice2D::cbegin() const { return data.cbegin(); }
Lattice2D::iterator Lattice2D::end() { return data.end(); }
Lattice2D::const_iterator Lattice2D::cend() const { return data.cend(); }


double Lattice2D::full_action(double lagrangian(Phi, Phi)) {
    int site, i;
    Phi forward_nphi_sum;
    vector<int> neighbors;
    double S = 0;
    // initial action
    for (int s=0; s<L*L; s++) {
        neighbors = neighbor_map[s];
        forward_nphi_sum  = data[neighbors[2]] + data[neighbors[3]];
        S += lagrangian(data[s], forward_nphi_sum);
    }

    return S;

}























