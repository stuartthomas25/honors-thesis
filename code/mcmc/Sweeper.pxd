cdef extern from "sweep.cpp":
    pass

cdef extern from "sweep.h" namespace "croutines":
    cdef cppclass Sweeper:
        Sweeper() except +
        Sweeper(int dim, double m02, double lam) except +
        int dim
        double redef_mass, quarter_lam
        int test(int i) except +
        int wrap(int c) except +
        void full_neighbors(int site, int neighbors[4]) except +
        double sweep(double data[], int sites[], int num_sites) except +
