cdef extern from "sweep.cpp" namespace "croutines":
    cdef cppclass Sweeper:
        Sweeper() except +
        int test() except +
