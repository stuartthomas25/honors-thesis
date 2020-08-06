# distutils: language = c++

'''
This file must be compiled to run. Call
    
    python setup.py build_ext --inplace

to compile this file. The functions can then be called through the `croutines` module. 
'''

import numpy as np
cimport numpy as np
cimport cython
from os import getpid   
from random import random

from libc.math cimport exp
from libc.stdlib cimport rand, RAND_MAX, srand
from libc.float cimport DBL_MAX
from libcpp.stack cimport stack
from libcpp.vector cimport vector
from libcpp.set cimport set

import sys

cdef extern from "<algorithm>" namespace "std" nogil:
    Iter find[Iter, T](Iter first, Iter last, const T& value) except +


def seed(int seed):
    '''Python wrapper to seed C rand()'''
    srand(seed)

cdef unsigned int wrap( int c, unsigned int dim):
    cdef int mod = c % <int>dim
    if mod < 0:
        mod += dim
    return mod
    
cdef void full_neighbors((int,int) site, (int, int) neighbors[4], int dim):    
    neighbors[0] = (wrap(site[0]+1,dim), site[1])
    neighbors[1] = (wrap(site[0]-1,dim), site[1])
    neighbors[2] = (site[0], wrap(site[1]+1,dim))
    neighbors[3] = (site[0], wrap(site[1]-1,dim))
    
# cdef void forward_neighbors((int,int) site, (int,int) neighbors[2], int dim):    
    # neighbors[0] = (site[0]+1, site[1])
    # neighbors[2] = (site[0], site[1]+1)

# cdef void backward_neighbors((int,int) site, (int,int) neighbors[2], int dim):    
    # neighbors[1] = (site[0]-1, site[1])
    # neighbors[3] = (site[0], site[1]-1)

cdef double randf():
    return <double>rand() / RAND_MAX

cdef int randint(int n):
    # Chop off all of the values that would cause skew...
    cdef int end = RAND_MAX / n # truncate skew
    end *= n

    cdef int r = rand()
    while r  >= end:
        r = rand()
    return r % n

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

    cdef double A
    cdef double r

    for i in range(num_sites):

        x = sites[i,0]
        y = sites[i,1]
        phi = data_view[x, y]
        
        backward_nphi_sum = data_view[(x - 1 + L) % L, y] + data_view[x, (y - 1 + L) % L] # add L here to ensure positive index
        forward_nphi_sum  = data_view[(x + 1) % L, y] + data_view[x, (y + 1) % L]
        
        dphi = crand_dist(randf())
        old_L = clagrangian( phi, forward_nphi_sum, redef_mass, quarter_lam)
        new_L = clagrangian( phi+dphi, forward_nphi_sum, redef_mass, quarter_lam)
        dS = (new_L - old_L) - backward_nphi_sum * dphi
        
        A = exp(-dS)
        r = randf()
        if dS < 0 or r <= A:
            dphis[i] = dphi
            tot_dS += dS

    return dphis_arr, tot_dS



### Generate Cluster


cdef double cPadd(double phi_a, double phi_b):
    return 1 - exp(-2*phi_a*phi_b) # Eq. 7.17

def Padd(phi_a, phi_b):
    return cPadd(phi_a, phi_b)

cdef int sign(double x):
    return (x > 0) - (x < 0)


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.cdivision(True)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef vector [(int,int)] generate_cluster(double [:,:] lat_data, (int,int) seed, int accept_all):
    cdef int phi_sign
    cdef double phi_a
    cdef double phi_b
    cdef (int,int) s
    cdef (int,int) n
    
    cdef int dim = lat_data.shape[0]
    
    cdef double Padd
    cdef stack [((int,int), double)] to_test
    cdef set [int] tested # for some reason, can't use a tuple for set
    cdef (int,int) neighbors[4]
    cdef vector [(int,int)] cluster
        
    phi_a = lat_data.__getitem__(seed)
    phi_sign = sign(phi_a) # get the sign of phi_a
    to_test.push((seed, DBL_MAX * phi_sign)) # site and phi value, using infinity to ensure Padd=1 for first addition
    
    while to_test.size()>0:
        s, phi_a = to_test.top()
        to_test.pop()
        if tested.find(s[0] + dim * s[1]) != tested.end(): # contains
            continue
        tested.insert( s[0] + dim * s[1] )
        
        phi_b = lat_data.__getitem__(s)
        if sign(phi_b) == phi_sign:
            Padd = cPadd(phi_a, phi_b)
            if accept_all or randf() < Padd: # check Padd==1 first to increase efficiency in SW algorithm
                cluster.push_back(s)
                full_neighbors(s, neighbors, dim)
                for n in neighbors:
                    if tested.find(n[0] + dim * n[1]) == tested.end(): #does not contain
                        to_test.push(( n, lat_data.__getitem__(s) ))

    return cluster

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.cdivision(True)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def wolff(lat):
    cdef int dim = lat.dim
    cdef (int, int) seed = (randint(dim), randint(dim))
    cdef vector [(int,int)] cluster
    cdef double[:,:] lat_view = lat.data
    cdef (int,int) neighbors[4]
    cdef (int,int) n

    cluster = generate_cluster(lat.data, seed, False)
    
    cdef double action = lat.action
    cdef double phi

    for c in cluster:
        phi = lat.data[c[0], c[1]]
        lat.data[c[0], c[1]] = -phi
        full_neighbors(c, neighbors, dim)
        for n in neighbors:
            action += 2 * lat_view[n[0], n[1]] * phi  #this derives as follows

    return lat.data, action



srand(int(random()*RAND_MAX)) # use python's random() to seed the terrible C rand() function.



