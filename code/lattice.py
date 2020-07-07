from functools import lru_cache
from itertools import product
from colorsys import hls_to_rgb
from random import randrange

from matplotlib import pyplot as plt
import numpy as np
from mpi4py import MPI

from keys import *

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()

REL_COORDS = ((-1,0),(0,-1),(1,0),(0,1))


class Lattice(object):

    def __init__(self, m, l, data=None, dim=None):
        self._construct( m, l, data, dim)


    def _construct(self, m, l, data, dim):
        self.m = m
        self.l = l

        self.__redefined_mass = 4 + 1/2 * self.m
        self.__quarter_lambda = 1/4 * self.l

        if data is not None:
            self.data = data
            self.dim = data.shape[0]
            self.size = data.size

        elif dim is not None:
            self.dim = dim
            self.size = dim**2
            if RANK==0:
                self.data = self.rand_dist(np.random.rand(self.dim, self.dim)) # Hot start
            else:
                self.data = np.empty((self.dim, self.dim))

        else:
            raise Exception("dim and data cannot both be None")

        self.__neighbors = lru_cache(maxsize=self.size) (self.__neighbors)
        self.sites = list(product(range(self.dim), repeat=2))
        # populate cache
        # for s in self.sites:
            # self.__neighbors(s)
        if RANK==0:
            self.calculate_action()
        else:
            self.action=0

        self.mpi_assignment = {}
        if RANK>0:
            global site_arrays
            for color in [WHITE, BLACK]:
                all_sites = list(self.checker_iter(color))
                split_sites = np.array_split(all_sites, SIZE-1)
                self.mpi_assignment[color] = split_sites[RANK-1]

    @staticmethod
    def rand_dist(r):
        return 3 * r - 1.5

    @staticmethod
    def static_lagrangian(phi, neighbor_phis, redef_mass, quarter_lam):
        return -phi * sum(neighbor_phis) + redef_mass * phi**2 + quarter_lam * phi**4

    def calculate_action(self):
        self.action = sum(self.lagrangian(c) for c in self.sites)

    def get_real_action(self):
        return sum(self.lagrangian(c) for c in self.sites)

    def __repr__(self):
        return repr(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getstate__(self):
        return self.m, self.l, self.data

    def __setstate__(self, state):
        self._construct(*state, None)

    def __copy__(self):
        return Lattice(self.m, self.l, np.copy(self.data))

    def __getitem__(self, k):
        return self.data[k]

    def __setitem__(self, k, v):
        self.data[k] = v

    def random(self):
        i = randrange(self.size)
        return divmod(i, self.dim)

    def checker_iter(self, color):
        offset = 1 if color==BLACK else 0
        for j in range(self.dim):
            for i in range((j+offset)%2, self.dim, 2):
                yield i,j

    def neighbors(self, coord):
        return self.__neighbors(tuple(coord))

    def __neighbors(self, coord):
        '''This function is cached in __init__'''
        d = self.dim
        x, y = coord
        return tuple( ((x+i) % self.dim, (y+j) % self.dim) for i,j in REL_COORDS )

    def lagrangian(self, coord):
        return self.static_lagrangian( self[coord], (self[c] for c in self.neighbors(coord)), self.__redefined_mass, self.__quarter_lambda )

    def show(self, figsize=(6,6), show=True):
        M = np.transpose(np.array(self.data))
        r = np.abs(M)
        arg = np.angle(M)

        # This code is generalized for complex fields
        h = (arg + np.pi)  / (2 * np.pi) + 0.5
        l = 1.0 - 1.0/(1.0 + r**1.0)
        s = 0.8

        c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
        c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
        c = c.swapaxes(0,2)
        plt.figure(figsize=figsize)
        im = plt.imshow(c)
        if show:
            plt.show()
        return im

    def magnetization(self):
        return np.sum(self.data) / self.size

    def susceptibility(self):
        m = self.magnetization()
        m2 = 1
        return (m2 - m**2)

    def binder_cumulant(self):
        return 0.
        phi_sq = np.sum(self.data**2)
        phi_qu = np.sum(self.data**4)
        return 1 - phi_qu / (3 * phi_sq**2)

