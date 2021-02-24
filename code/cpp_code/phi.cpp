#include "phi.h"

using namespace std;

Phi::Phi() {}; 

Phi::Phi(array<double, N> phi){
    this->phi = phi;
};

void Phi::init_as_zero(){
    phi = {0,0,0};
    //for (int i=0; i<N; i++) {
        //phi[i] = 0;
    //}
}

double Phi::norm_sq() const {
    double cumsum = 0;
    for (int i=0; i<N; i++) {
        cumsum += phi[i] * phi[i];
    }
    return cumsum;

}

Phi& Phi::operator+=(const Phi& other) {
    for (int i=0; i<N; i++) {
        phi[i] += other[i];
    }
    return *this;
}

Phi& Phi::operator*=(const double & a) {
    for (int i=0; i<N; i++) {
        phi[i] *= a;
    }
    return *this;
}

Phi& Phi::operator/=(const double & a) {
    for (int i=0; i<N; i++) {
        phi[i] /= a;
    }
    return *this;
}

Phi Phi::operator+ (const Phi & aPhi) const {
    Phi new_phi(phi);
    //for (int i=0; i<N; i++) {
        //new_phi[i] = phi[i] + aPhi[i];
    //}
    new_phi += aPhi;
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

Phi Phi::operator& (const Phi & other) const {
    Phi new_phi;
    new_phi[0] = phi[1] * other[2] - phi[2] * other[1];
    new_phi[1] = phi[2] * other[0] - phi[0] * other[2];
    new_phi[2] = phi[0] * other[1] - phi[1] * other[0];
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
    Phi new_phi(phi);
    new_phi *= a;
    return new_phi;
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

Phi operator*(double a, const Phi& b)
{
    return b*a;
}


