from random import random
from math import exp
from keys import WOLFF, SWENDSEN_WANG

import numpy as np
from mpi4py import MPI

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()

if 'line_profiler' not in globals():
    def profile(f):
        return f


class RandomWalk(object):
    def __init__(self, lattice):
        self.lat = lattice
        self.mpi_chunksize = lattice.size ** 2 // 2 // SIZE // 2

    @profile
    def checker_met(self, site, color):
        neighbor_phis = tuple(self.lat.data[n] for n in self.lat.neighbors(site)[:2])
        old_action = self.lat.action
        old_L = self.lat.lat_lagrangian(site)
        dphi = self.lat.rand_dist(random())
        self.lat.data[site] += dphi
        new_L = self.lat.lat_lagrangian(site)
        dS = new_L - old_L
        for nphi in neighbor_phis:
            dS -= nphi * dphi

        if dS < 0 or random() <= exp(-dS): # EDIT
            return dphi, dS
        else:
            self.lat.data[site] -= dphi
            self.lat.action = old_action
            return None


    @profile
    def half_checkerboard(self, color):
        COMM.Bcast(self.lat.data, root=0)
        self.lat.action = COMM.bcast(self.lat.action, root=0)

        # One node compiles changes from other nodes
        if RANK==0:
            running_nodes = SIZE - 1
            while running_nodes > 0:
                data = COMM.recv()
                xs = data['xs']
                ys = data['ys']
                dphis = data['dphis']
                dS = data['dS']

                self.lat.data[xs, ys] += dphis
                self.lat.action += dS
                if data['done']:
                    running_nodes -= 1

        else:
            size = 0
            xs = []
            ys = []
            dphis = []
            tot_dS = 0
            for site in self.lat.mpi_assignment[color]:
                res = self.checker_met(tuple(site), color)
                if res is not None:
                    dphi, dS = res
                    xs.append(site[0])
                    ys.append(site[1])
                    dphis.append(dphi)
                    tot_dS += dS
                    size += 1

                    if size % self.mpi_chunksize == 0:
                        COMM.send({'xs':xs, 'ys':ys, 'dphis':dphis, 'dS':dS, 'done':False}, dest=0)
                        xs = []
                        ys = []
                        dphis = []
                        tot_dS = 0


            COMM.send({'xs':xs, 'ys':ys, 'dphis':dphis, 'dS':tot_dS, 'done':True}, dest=0)


    def checkerboard(self):
        if SIZE<2:
            raise Exception("Please run with more than one core")
        self.half_checkerboard(False)
        self.half_checkerboard(True)

    def wolff(self):
        seed = self.lat.random()
        cluster = self.generate_cluster(seed, False)

        if len(cluster) <= self.lat.size // 2:
            for c in cluster:
                self.lat.data[c] *= -1
        else:
            for c in self.lat.sites:
                if not tuple(c) in cluster:
                    self.lat.data[c] *= -1

        self.lat.calculate_action()

    @staticmethod
    def Padd(phi_a, phi_b):
        Padd = 1 - exp(-2*phi_a*phi_b) # Eq. 7.17
        assert Padd>0
        return Padd

    def generate_cluster(self, seed, accept_all):
        sign = np.sign(self.lat[seed])
        to_test = {(seed, np.inf*sign)} # site and phi value, using infinity to ensure Padd=1 for first addition
        tested = set()
        cluster = set()
        while len(to_test)>0:
            s, phi_a = to_test.pop()
            tested.add(s)
            if np.sign(self.lat[s]) == sign:
                Padd = self.Padd(phi_a, self.lat[s])
                if Padd==1 or random()<Padd: # check Padd==1 first to increase efficiency in SW algorithm
                    cluster.add(s)
                    to_test |= set((n, self.lat[s]) for n in self.lat.neighbors(s) if n not in tested)

        return cluster

    def run(self,
            sweeps,
            cluster_method=None,
            cluster_rate=5,
            record_rate=5,
            thermalization=0,
            recorder=None
            ):

        if cluster_method==WOLFF:
            cluster = self.wolff
        elif cluster_method==SWENDSEN_WANG:
            cluster = self.swendsen_wang
        else:
            cluster = lambda: None

        if thermalization==0 and recorder is not None:
            recorder.record(self.lat)

        for i in range(sweeps):
            self.checkerboard()

            if recorder is not None and i % record_rate == 0 and i > thermalization:
                recorder.record(self.lat)
            if i % cluster_rate ==0 and RANK==0:
                cluster()

               #  if recorder is not None and i % record_rate == 0 and i > thermalization:
                    # recorder.save(self.lat)


