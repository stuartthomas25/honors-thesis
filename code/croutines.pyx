'''
This file must be compiled to run. Call
    
    python setup.py build_ext --inplace

to compile this file. The functions can then be called through the `croutines` module. 
'''


import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp
from libc.stdlib cimport rand, RAND_MAX, srand
from os import getpid   
from random import random

def seed(int seed):
    '''Python wrapper to seed C rand()'''
    srand(seed)
    
cdef double clagrangian(double phi, double nphi_sum, double redef_mass, double quarter_lam):
    return -phi * nphi_sum + redef_mass * phi**2 + quarter_lam * phi**4
    
# cdef double crand_dist(double r):
    # return 3 * r - 1.5
cdef crand_dist(double r):
    return 3 * r - 1.5

def lagrangian(double phi, double nphi_sum, double redef_mass, double quarter_lam):
    return clagrangian(phi, nphi_sum, redef_mass, quarter_lam)

def rand_dist(double r):
    return crand_dist(r)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.cdivision(True)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def sweep(lat, int [:,:] sites):
    cdef int L = lat.dim
    cdef double redef_mass = 2 + 0.5 * lat.m02
    cdef double quarter_lam = 0.25 * lat.lam
    cdef double [:, :] data_view = lat.data
    cdef int num_sites = len(sites)
    
    dphis_arr = np.zeros(num_sites, dtype='d')
    cdef double [:] dphis = dphis_arr
    cdef double tot_dS = 0.
    cdef double dS, new_L, old_L, dphi, phi, backward_nphi_sum, forward_nphi_sum
    cdef int i, x, y

    for i in range(num_sites):

        x = sites[i,0]
        y = sites[i,1]
        phi = data_view[x, y]

        backward_nphi_sum = data_view[(x - 1 + L) % L, y] + data_view[x, (y - 1 + L) % L] # add L here to ensure positive index
        forward_nphi_sum  = data_view[(x + 1) % L, y] + data_view[x, (y + 1) % L]
        
        dphi = crand_dist(<double>rand() / RAND_MAX)
        old_L = clagrangian( phi, forward_nphi_sum, redef_mass, quarter_lam)
        new_L = clagrangian( phi+dphi, forward_nphi_sum, redef_mass, quarter_lam)
        dS = (new_L - old_L) - backward_nphi_sum * dphi

        if dS < 0 or <double>rand() / RAND_MAX <= exp(-dS):
            dphis[i] = dphi
            tot_dS += dS

    return dphis_arr, tot_dS


srand(int(random()*RAND_MAX)) # use python's random() to seed the terrible C rand() function.
