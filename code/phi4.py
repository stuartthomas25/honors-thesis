from random import randint, choices, random, seed
from colorsys import hls_to_rgb
from copy import copy, deepcopy
import numpy as np
import time
from functools import lru_cache
from matplotlib import pyplot as plt
import sys
import imageio

class Lattice(object):

    rel_coords = [(-1,0),(0,-1),(1,0),(0,1)]

    def __init__(self, m, l, data=None, dim=None):
        self._construct( m, l, data, dim)

    def _construct(self, m, l, data, dim):
        self.m = m
        self.l = l

        if data:
            self.data = data
            self.dim = int(np.sqrt(len(data)))
            self.size = len(data)

        elif dim:
            self.dim = dim
            self.size = dim**2
            self.data = [self._rand_init_val() for _ in range(self.size)]

        else:
            raise Exception("dim and data cannot both be None")

        self.sites = list(range(self.size))
        self.energy = 0
        for coord in range(self.size):
            self.energy += self.site_hamiltonian(coord)

        self.neighbors = lru_cache(maxsize=self.size) (self.neighbors)

    def _rand_init_val(self):
        # random value in [-1.5, 1.5)
        return random() * 3 - 1.5 # EDIT

    def __repr__(self):
        return repr(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, val):
        self.data[i] = val
        self.update_energy(i)

    def __getstate__(self):
        return self.m, self.l, self.data

    def __setstate__(self, state):
        self._construct(*state, None)

    def random(self):
        return randint(0,self.size-1)

    def coords(self):
        return range(self.size)

    def neighbors(self, coord):
        d = self.dim
        x = coord % d
        y = coord // d
        ns = []
        for i,j in Lattice.rel_coords:
            ns.append( (x+i) % self.dim + d*((y+j) % self.dim) )
        return tuple(ns)

    def site_hamiltonian(self, site):
        e = 0
        phi = self[site]
        for neighbor in self.neighbors(site):
            e -= phi * self[neighbor]

        e += 0.5 * self.m**2 * phi**2 + 0.25 * self.l * phi**4 # based on 7.16, QUESTION

        return e

    def update_energy(self, *sites):
        for site in sites:
            self.energy += 4 * self.site_hamiltonian(site)

    def show(self, figsize=(6,6), show=True):
        M = np.transpose(np.resize( np.array(self.data), (self.dim, self.dim) ))
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
        return sum(self) / self.size

    def susceptibility(self):
        m = self.magnetization()
        m2 = 1
        return (m2 - m**2)

    def binder_cumulant(self):
        phi_sq = sum(phi**2 for phi in self)
        phi_qu = sum(phi**4 for phi in self)
        return 1 - phi_qu / (3 * phi_sq**2)

class RandomWalk(object):
    def __init__(self, lat):
        self.lattice = lat

    def flip(self, i):
        self.lattice[i] =-self.lattice[i]

    def change(self, i):
        self.flip(i)
        self.lattice[i] += 3 * random() - 1.5

    def metropolis(self):
        '''Returns True if accepted'''
        e = self.lattice.energy
        site = self.lattice.random()
        old_val = self.lattice[site]
        self.change(site)
        new_e = self.lattice.energy

        if new_e > e:
            A = np.exp( e - new_e )
            accept = random() < A
            if not accept:
                self.lattice[site] = old_val
                return False
        return True

    def wolff(self):
        seed = self.lattice.random()
        cluster = self._generate_cluster(seed, False)

        if len(cluster) <= self.lattice.size // 2:
            for i in cluster:
                self.flip(i)
        else:
            for i in self.lattice.coords():
                if not i in cluster:
                    self.flip(i)

    def swendsen_wang(self):
        sites_left = set(self.lattice.sites)
        while len(sites_left) > 0:
            seed = sites_left.pop()
            cluster = self._generate_cluster(seed, True)
            sites_left -= cluster
            if random() < 0.5:
                for i in cluster:
                    self.flip(i)

    def _Padd(self, site_out, site_in):
        return 1 - np.exp(-2*self.lattice[site_in]*self.lattice[site_out])


    def _generate_cluster(self, seed, accept_all):
        sign = np.sign(self.lattice[seed])
        to_test = {(seed, 1)} # site and Padd
        tested = set()
        cluster = set()
        while len(to_test)>0:
            s, Padd = to_test.pop()
            tested.add(s)
            if np.sign(self.lattice[s]) == sign:
                if Padd==1 or random()<Padd: # check Padd==1 first to increase efficiency in SW algorithm
                    cluster.add(s)
                    to_test |= set((n, self._Padd(n, s)) for n in self.lattice.neighbors(s) if n not in tested)

        return cluster




class GifProducer(object):
    def __init__(self):
        self.frames = []

    def save_lat(self, lat):
        im = lat.show(show=False)
        self.frames.append(im.make_image("AGG")[0])
        plt.close()

    def save(self, fn, fps=3):
        imageio.mimwrite(fn, self.frames, fps=fps)


if __name__=="__main__":
    seed(a=14231)
    L = 256
    m = 1
    lam = 1
    wolff = False
    l = Lattice(dim=L, m=m, l=lam)
    rw = RandomWalk(l)
    wolff_rate = 5
    stop_acceptance = 0.001

    counter = 0
    magnetizations = []
    energies = []
    susceptibilities = []
    binder_cums = []
    states = []
    sweep_count = 0
    record = True
    record_count = 0
    record_rate = l.size
    thermalization = 200
    iter_index = []
    gp = GifProducer()
    def record_state(lat):
        iter_index.append(sweep_count)
        susceptibilities.append(l.susceptibility())
        binder_cums.append(l.binder_cumulant())
        magnetizations.append(abs(l.magnetization()))
        energies.append(l.energy/l.size)
        gp.save_lat(l)

    record_state(l)
    timeout = 100

    start = time.time()
    running = True
    while running and sweep_count<timeout:
        for _ in range(wolff_rate):
            percent_accepted = 0
            for _ in range(l.size):
                percent_accepted += rw.metropolis() / l.size
                record_count += 1
                if record and record_count > thermalization and record_count % record_rate == 0:
                    states.append(deepcopy(l))
            sweep_count += 1
            record_state(l)

            if percent_accepted < stop_acceptance:
                running = False
        if wolff:
            rw.wolff()
        else:
            rw.swendsen_wang()

        record_state(l)
    exec_time = time.time() - start
    print(f"Done in {exec_time}")

    gp.save("test.gif")

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(18,12), sharex = True)
    ax1.plot(iter_index, magnetizations)
    ax2.plot(iter_index, energies)
    ax3.plot(iter_index, susceptibilities)
    ax4.plot(iter_index, binder_cums)

    ax1.set_ylabel("Magnetization")
    ax2.set_ylabel("Total Energy")
    ax3.set_ylabel("Susceptibility")
    ax4.set_ylabel("Binder Cumulant")
    ax4.set_xlabel("Sweep")
    # ax3.set_ylim((-0.05,1.05))
    ax1.set_title(f"Monte Carlo Simulation of Ising Model using Metropolis and {'Wolff' if wolff else 'Swenson-Wang'} Algorithms, $L={L}$, $\\lambda={lam}$, $\\mu_0^2={m}$, $t={exec_time:.1f}s$")
    plt.savefig('plots/phi_248_sw.png')
    plt.show()


