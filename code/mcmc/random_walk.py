from random import random
from math import exp
from .keys import WOLFF, SWENDSEN_WANG, WHITE, BLACK
from .croutines import sweep, generate_cluster, wolff
from time import time
import sys
import cProfile, pstats



import numpy as np
from mpi4py import MPI


CHECK_ACTION = False

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()

# if 'profile' not in globals():
    # profile = lambda x : x
# else:
    # print("profiling...")



class RandomWalk(object):
    def __init__(self, lattice):
        self.lat = lattice
        self.mpi_chunksize = lattice.size ** 2 // 2 // SIZE // 2

        # assign sites for MPI parallelizing
        self.mpi_assignments = [{} for r in range(SIZE)]

        for color in [WHITE, BLACK]:
            all_sites = list(self.lat.checker_iter(color))
            if len(all_sites) % SIZE != 0:
                raise Exception("L should be divisible by MPI group size")

            split_sites = np.array_split(all_sites, SIZE)
            for r in range(SIZE):
                self.mpi_assignments[r][color] = split_sites[r]


    def half_checkerboard(self, color):
        data_bk = np.copy(self.lat.data)
        COMM.Bcast(self.lat.data, root=0)
        self.lat.action = COMM.bcast(self.lat.action, root=0)

        sites = self.mpi_assignments[RANK][color]

        start = time()
        dphis, dS = sweep(self.lat, sites.astype('i'))

        if RANK==0:
            recv_dphis = np.empty([SIZE, dphis.size], dtype=dphis.dtype)
        else:
            recv_dphis = None

        COMM.Gather(dphis, recv_dphis, root=0)
        dS = COMM.gather(dS, root=0)

        if RANK==0:
            self.lat.data = data_bk
            for r in range(SIZE):
                sites = self.mpi_assignments[r][color]
                self.lat.data[sites[:,0], sites[:,1]] += recv_dphis[r]

            self.lat.action += sum(dS)

    def checkerboard(self):
        self.half_checkerboard(False)
        self.half_checkerboard(True)

    def flip(self, c):
        phi = self.lat.data[c]
        self.lat.data[c] = -phi
        self.lat.action += 2 * phi * sum(self.lat.data[n] for n in self.lat.half_neighbors[c])


    def wolff(self):
        seed = self.lat.random()
        # print(self.lat.data.shape, seed)
        cluster = generate_cluster(self.lat.data, seed, False)

        if len(cluster) <= self.lat.size // 2:
            for c in cluster:
                self.flip(c)
        else:
            for c in self.lat.sites:
                if not tuple(c) in cluster:
                    self.flip(c)

  #   @staticmethod
    # def Padd(phi_a, phi_b):
        # Padd = 1 - exp(-2*phi_a*phi_b) # Eq. 7.17
        # assert Padd>0
        # return Padd

    # def generate_cluster(self, seed, accept_all):
        # sign = np.sign(self.lat[seed])
        # to_test = {(seed, np.inf*sign)} # site and phi value, using infinity to ensure Padd=1 for first addition
        # tested = set()
        # cluster = set()
        # while len(to_test)>0:
            # s, phi_a = to_test.pop()
            # tested.add(s)
            # if np.sign(self.lat[s]) == sign:
                # Padd = self.Padd(phi_a, self.lat[s])
                # if Padd==1 or random()<Padd: # check Padd==1 first to increase efficiency in SW algorithm
                    # cluster.add(s)
                    # to_test |= set((n, self.lat[s]) for n in self.lat.full_neighbors[s] if n not in tested)

        # return cluster

    def run(self,
            sweeps,
            cluster_method=None,
            cluster_rate=5,
            hooks=None,
            recorders=[],
            progress=False,
            record_on_cluster=False
            ):

        progress_char = '#'

        start = time()

        if cluster_method==WOLFF:
            cluster = wolff
        else:
            cluster = None
        if hooks is None:
            hooks = [lambda x : x for r in recorders]
        if len(hooks) != len(recorders):
            raise Exception("Hooks and recorders don't match!")

        if progress and RANK==0:
            progressbar_size = 140
            progressbar_rate = int(sweeps/progressbar_size)
            print('['+'-'*progressbar_size+']', end='\r', flush=True)

        for i in range(sweeps):
            self.checkerboard()
            if RANK==0:
                if i % cluster_rate ==0 and cluster is not None:
                    self.lat.data, self.lat.action = cluster(self.lat)

                for j,r in enumerate(recorders):
                    if i % r.rate == 0 and i >= r.thermalization:
                        r.record( hooks[j](self.lat) )

                # for j,r in enumerate(recorders):
                    # if i % r.rate == 0 and i >= r.thermalization:
                        # evolved_lat = hooks[j](self.lat)
                        # print("S", np.sum(evolved_lat.data) / evolved_lat.size)
                        # r.record(evolved_lat)
                        # print("V", r.values['phi'][-1])
                        # hooks[j](self.lat).show()

                if progress and i % progressbar_rate == 0:
                    p = int(i/sweeps*progressbar_size)
                    print('['+progress_char*p+'-'*(progressbar_size-p)+']', end='\r', flush=True)

        if progress and RANK==0:
            print('['+progress_char*progressbar_size+']',flush=True)


        for r in recorders:
            r.execution_time = time() - start
            r.finalize_values()



    def check_action(self, tol=1e5):
        assert np.isclose(self.lat.action, self.lat.calculate_action())

