from time import time, sleep
import sys
import os
from getopt import getopt
from mpi4py import MPI

from lattice import Lattice, show
from random_walk import RandomWalk
from keys import *
from recorder import Recorder
import numpy as np
from parallel_utils import COMM, RANK, SIZE


if 'line_profiler' not in globals():
    def profile(f):
        return f

MPI_RES_CHUNKSIZE = 0

def main():
    opts, args = getopt(sys.argv[1:],'pL:s:S')
    dopts = dict(opts)

    if '-s' in dopts:
        seed = int(dopts['-s'])
        seed(a=seed)
        np.random.seed(seed)

    L = int(dopts.get('-L', 64))
    m02 = -0.68
    lam = 0.5
    sweeps = 100
    cluster_method = None



    l = Lattice(dim=L, m02=m02, l=lam)



# Test gradient flow
#     recorder = Recorder()

    # def laplacian(a):
        # return np.gradient(np.gradient(a, axis=0), axis=0) + np.gradient(np.gradient(a, axis=1), axis=1)

    # last_rho = None
    # for tau in np.linspace(0.,0.005,100):
        # rho = l.rho(tau)
        # recorder.save(rho)
        # if last_rho is not None:
            # dt_data = rho.data - last_rho.data
            # grad_data = laplacian(rho.data)
        # last_rho = rho
    # recorder.plot(title=f"Effect of gradient flow on observables: $L={L}$, $\\lambda={lam}$, $\\mu_0^2={m02}$", fname="gradient_flow", xlabel=r'$\tau$ (flow time)', show=True)

    # recorder.save_gif("plots/gradient_flow.gif", fps=30)

    # quit() # just testing the gradient flow, so stop here

    #set the site assignments (factor of 2 for checkerboard)


    cluster_rate = 1
    record_rate = 1
    thermalization = 1

    recorder = Recorder()
    recorder.save(l)

    rw = RandomWalk(l)

    start = time()

    rw.run(sweeps, cluster_method, cluster_rate, record_rate, thermalization, recorder)

    if RANK==0:
        exec_time = time() - start
        print(f"{exec_time}")

        recorder.save_gif("plots/lattice_visualization.gif")

        recorder.plot(title=f"Monte Carlo Simulation of $\phi^4$ Model using Metropolis and Wolff Algorithms, $L={L}$, $\\lambda={lam}$, $\\mu_0^2={m02}$, $t={exec_time:.1f}s$", fname="parallel.png", show=('-S' in dopts))

if __name__=="__main__":

    main()
