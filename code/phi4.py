# from random import randint, choices, random, seed, randrange
# from copy import copy, deepcopy
# import numpy as np
from time import time, sleep
# from functools import lru_cache
# from itertools import product
import sys
# import imageio
# import timeit
# from math import exp, sqrt
# import inspect
# import tempfile
# import multiprocessing as mp
# import pickle
import os
from getopt import getopt
# import tables
# from dataclasses import dataclass
# from typing import Callable
# import dill
# from matplotlib import pyplot as plt
from mpi4py import MPI
# from datetime import datetime

from lattice import Lattice
from random_walk import RandomWalk
from keys import *
from recorder import Recorder


COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


if 'line_profiler' not in globals():
    def profile(f):
        return f

MPI_RES_CHUNKSIZE = 0

def main():
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

    opts, args = getopt(sys.argv[1:],'pL:s:S')
    dopts = dict(opts)

    if '-s' in dopts:
        seed = int(dopts['-s'])
        seed(a=seed)
        np.random.seed(seed)

    L = int(dopts.get('-L', 8))
    m = 1
    lam = 1
    iterations = 1
    cluster_method = WOLFF



    l = Lattice(dim=L, m=m, l=lam)

    #set the site assignments (factor of 2 for checkerboard)


    wolff_rate = 1

    recorder = Recorder()
    recorder.save(l)

    rw = RandomWalk(l)

    start = time()
    for _ in range(iterations):
        for _ in range(wolff_rate):
            rw.checkerboard()
            recorder.save(l)

        if RANK==0:
            if cluster_method==WOLFF:
                rw.wolff()
            elif cluster_method==SWENSDEN_WANG:
                rw.swendsen_wang()

        recorder.save(l)

    if RANK==0:
        exec_time = time() - start
        print(f"{exec_time}")

        recorder.save_gif("plots/lattice_visualization.gif")

        recorder.plot(L, m, lam, exec_time, '-S' in dopts)


if __name__=="__main__":
    main()
