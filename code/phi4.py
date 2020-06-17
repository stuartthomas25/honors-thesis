from random import randint, choices, random
from colorsys import hls_to_rgb
from copy import copy as copy
import numpy as np
import time
from functools import lru_cache
from matplotlib import pyplot as plt
import sys

class Lattice(object):

    rel_coords = [(-1,0),(0,-1),(1,0),(0,1)]

    def __init__(self, dim, beta, J):
        self.beta = beta
        self.dim = dim
        self.size = dim**2
        self.data = choices([True,False], k=self.size)
        self.J = J
        self.energy = 0
        for coord in range(self.size):
            self.energy += self.site_hamiltonian(coord)

    def __repr__(self):
        return repr(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, val):
        self.data[i] = val
        self.update_energy(i)

    def random(self):
        return randint(0,self.size-1)

    def coords(self):
        return range(self.size)

    def neighbors(self, coord):
        d = self.dim
        x = coord % d
        y = coord // d

        for i,j in Lattice.rel_coords:
            if 0 <= x+i < self.dim and 0 <= y+j < self.dim:
                yield (x+i) + d*(y+j)

    def site_hamiltonian(self, site):
        e = 0
        for neighbor in self.neighbors(site):
                if self[site] ^ self[neighbor]:
                    e +=  1
                else:
                    e += -1

        return self.J * e

    def update_energy(self, *sites):
        for site in sites:
            self.energy += 4 * self.site_hamiltonian(site)

    def show(self, figsize=(6,6)):
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
        plt.imshow(c)
        plt.show()

    def magnetization(self):
        return sum([1 if phi==True else -1 for phi in self]) / self.size

    def susceptibility(self):
        m = self.magnetization()
        m2 = 1
        return self.size * self.beta * (m2 - m**2)

class RandomWalk(object):
    def __init__(self, lat):
        self.lattice = lat
        self.beta = lat.beta
        self.J = lat.J
        self.Padd = 1 - np.exp(-2 * self.beta * self.J)

    def flip(self, i):
        self.lattice[i] = not self.lattice[i]

    def metropolis(self):
        '''Returns True if accepted'''
        e = self.lattice.energy
        site = self.lattice.random()
        self.flip(site)
        new_e = self.lattice.energy

        if new_e > e:
            A = np.exp( -self.beta * (new_e - e) )
            accept = random() < A
            if not accept:
                self.flip(site)
                return False
        return True

    def wolff(self):
        seed = self.lattice.random()
        old_recursion_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(self.lattice.size)
        cluster = set([seed])
        self._generate_cluster(cluster, seed, self.lattice[seed], 0)
        if len(cluster) >= self.lattice.size // 2:
            for i in cluster:
                self.flip(i)
        else:
            for i in self.lattice.coords():
                if not i in cluster:
                    self.flip(i)

        sys.setrecursionlimit(old_recursion_limit)


    def _generate_cluster(self, cluster, seed, value, i):
        neighbors = set(self.lattice.neighbors(seed)) - cluster
        cluster.update(neighbors)
        for neighbor in neighbors:
            if self.lattice[neighbor]==value:
                if random() < self.Padd:
                    self._generate_cluster(cluster, neighbor, value, i+1)


if __name__=="__main__":
    L = 128
    beta = 100
    J = 1
    l = Lattice(L, J, beta)
    rw = RandomWalk(l)
    wolff_rate = 5
    stop_acceptance = 0.001

    counter = 0
    magnetizations = []
    energies = []
    susceptibilities = []
    sweep_count = 0
    iter_index = []
    def record_state(lat):
        iter_index.append(sweep_count)
        susceptibilities.append(l.susceptibility())
        magnetizations.append(abs(l.magnetization()))
        energies.append(l.energy)

    record_state(l)

    start = time.time()
    running = True
    while running:
        for _ in range(wolff_rate):
            percent_accepted = 0
            for _ in range(l.size):
                percent_accepted += rw.metropolis() / l.size
            sweep_count += 1
            record_state(l)

            if percent_accepted < stop_acceptance:
                running = False

        rw.wolff()
        record_state(l)

    print(f"Done in {time.time() - start}")

    fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = True)
    ax1.plot(iter_index, magnetizations)
    ax2.plot(iter_index, energies)
    ax3.plot(iter_index, susceptibilities)

    ax1.set_ylabel(r"Avg. Magnetization $\|\langle m \rangle\|$")
    ax2.set_ylabel("Total Energy $E$")
    ax3.set_ylabel("Susceptibility $\chi$")
    ax3.set_xlabel("Sweep")
    ax1.set_title(f"Monte Carlo Simulation of Ising Model using Metropolis and Wolff Algorithms, $L={L}$, $\\beta={beta}$, $J={J}$")

    plt.show()


