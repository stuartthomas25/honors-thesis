from functools import lru_cache
from itertools import product
from colorsys import hls_to_rgb
from random import randrange
import abc
from .croutines import lagrangian, rand_dist
from copy import copy

import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

from matplotlib import pyplot as plt
import numpy as np
from mpi4py import MPI

from .keys import *

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()

REL_COORDS = ((1,0),(0,1),(-1,0),(0,-1))

def show(data, figsize=(6,6), show=True):
    plt.figure(figsize=figsize)
    im = plt.imshow(c)
    if show:
        plt.show()
    return im

class Lattice(object, metaclass=abc.ABCMeta):

    def __init__(self, data=None, dim=None):
        self._construct(data, dim)


    def _construct(self, data, dim):
        if data is not None:
            self.data = data
            self.dim = data.shape[0]
            self.size = data.size


        elif dim is not None:
            self.dim = dim
            self.size = dim**2
            if RANK==0:
                self.data = np.vectorize(rand_dist) (np.random.rand(self.dim, self.dim)) # Hot start
            else:
                self.data = np.empty((self.dim, self.dim))

        else:
            raise Exception("dim and data cannot both be None")

        self.sites = list(product(range(self.dim), repeat=2))
        self.full_neighbors = {s:self.neighbors(s) for s in self.sites}
        self.half_neighbors = {s:self.neighbors(s)[:2] for s in self.sites}

        if RANK==0:
            self.action = self.calculate_action()
        else:
            self.action = 0

    def calculate_action(self):
        action = 0
        for c in self.sites:
            action += self.lat_lagrangian(c)
        return action

    def __repr__(self):
        return repr(self.data)

    def __iter__(self):
        return iter(self.data)

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

    @abc.abstractmethod
    def lat_lagrangian(self, coord):
        pass

    def neighbors(self, coord):
        d = self.dim
        x, y = coord
        return tuple( ((x+i) % self.dim, (y+j) % self.dim) for i,j in REL_COORDS )


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
        m2 = np.sum(self.data**2) / self.size
        return (m2 - m**2)

    def binder_cumulant(self):
        phi_sq = np.sum(self.data**2) / self.size
        phi_qu = np.sum(self.data**4) / self.size
        return 1 - phi_qu / (3 * phi_sq**2)



class Phi4Lattice(Lattice):

    def __init__(self, m02, lam, dim):
        self.m02 = m02
        self.lam = lam

        self.redef_mass = 2 + 0.5 * m02
        self.quarter_lam = 0.25 * lam

        super().__init__(dim=dim)

    def lat_lagrangian(self, coord):
        return lagrangian( self[coord], sum(self[c] for c in self.neighbors(coord)[:2]), self.redef_mass, self.quarter_lam)

class GradientFlow(object):
    def flow_evolution(self, lat, tau, action=False):
        fft = np.fft.fft2(lat.data)
        ps = np.concatenate([np.arange(lat.dim//2), -np.arange(lat.dim//2,0,-1)])
        px, py = np.meshgrid(ps, ps)
        new_lat = copy(lat) # make shallow copy so as not to copy the entire field
        new_lat.data = np.exp(-tau*(px**2 + py**2)) * fft
        if action:
            new_lat.calculate_action()
        return new_lat
