from random import randint, choices
from colorsys import hls_to_rgb
from copy import copy as copy
import numpy as np
import time
from functools import lru_cache
from matplotlib import pyplot as plt

class Lattice(object):

    rel_coords = [(-1,0),(0,-1),(1,0),(0,1)]

    def __init__(self, dim):
        self.dim = dim
        self.size = dim**2
        self.data = choices([True,False], k=self.size)

    def __repr__(self):
        return repr(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, val):
        self.data[i] = val

    def random(self):
        return randint(0,self.size-1)

    def coords(self):
        return range(self.size)

    def neighbors(self, coord):
        d = self.dim
        x = coord % d
        y = coord // d

        for i,j in Lattice.rel_coords:
            if 0 < (x+i) + d*(y+j) < self.size:
                yield (x+i) + d*(y+j)

    def site_hamiltonian(self, site):
        for neighbor in self.neighbors(site):
                if self[site] ^ self[neighbor]:
                    return 1
                else:
                    return -1

    def hamiltonian(self):
        energy = 0
        J = 1
        for coord in range(self.size):
            energy += self.site_hamiltonian(coord)
        return energy

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


class RandomWalk(object):
    def __init__(self, lat, beta):
        self.lattice = lat
        self.beta = beta

    def random_change(self):
        i = self.lattice.random()
        self.flip(i)
        return i

    def flip(self, i):
        self.lattice[i] = not self.lattice[i]

    def metropolis(self):
        '''Returns True if accepted'''
        e = self.lattice.hamiltonian()
        old_lat = copy(self.lattice.data)
        site = self.random_change()
        new_e = self.lattice.hamiltonian()

        if new_e > e:
            A = np.exp( -self.beta * (new_e - e) )
            accept = choices([True,False], cum_weights=[A,1])[0]
            if not accept:
                self.lattice.data = old_lat
                return False
        return True

    def wolff(self):
        seed = self.lattice.random()
        cluster = set([seed])
        self.generate_cluster(cluster, seed, self.lattice[seed])
        if len(cluster) >= self.lattice.size // 2:
            for i in cluster:
                self.flip(i)
        else:
            for i in self.lattice.coords():
                if not i in cluster:
                    self.flip(i)

    def generate_cluster(self, cluster, seed, value):
        for neighbor in self.lattice.neighbors(seed):
            if not neighbor in cluster and self.lattice[neighbor]==value:
                cluster.add(neighbor)
                self.generate_cluster(cluster, neighbor, value)

if __name__=="__main__":
    L = 32
    beta = 1e8
    l = Lattice(L)
    rw = RandomWalk(l, beta)
    wolff_rate = l.size * 2

    rejections = 0
    counter = 0
    magnetizations = []
    energies = []

    while rejections < 100:
        if counter % 500 == 0:
            #l.show()
            pass
        accept = rw.metropolis()
        if accept:
            rejections = 0
        else:
            rejections += 1

        if counter % wolff_rate == 0:
            rw.wolff()


        magnetizations.append(l.magnetization())
        energies.append(l.hamiltonian())

        counter += 1


    print("Done")
    l.show()
    iter_index = list(range(len(magnetizations)))
    wolff_positions = list(range(0, len(magnetizations), wolff_rate))

    fig, (ax1, ax2) = plt.subplots(2,1, sharex = True)
    for pos in wolff_positions:
        ax1.axvline(pos,c='r',ls='--')
        ax2.axvline(pos,c='r',ls='--')
    ax1.plot(iter_index, magnetizations)
    ax2.plot(iter_index, energies)

    ax1.set_ylabel("Avg. Magnetization")
    ax2.set_ylabel("Total Energy")
    ax2.set_xlabel("Metropolis iteration")
    ax1.set_title(f"Monte Carlo Simulation of Ising Model using Metropolis and Wolff Algorithms, $L={L}$, $\\beta={beta}$")

    plt.show()

    print(l.hamiltonian())

