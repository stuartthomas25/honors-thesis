from parallel_utils import RANK
from keys import WOLFF, SWENDSEN_WANG
from lattice import Phi4Lattice
from random_walk import RandomWalk
from matplotlib import pyplot as plt
from recorder import Recorder
import numpy as np
from time import time
from copy import deepcopy

def main():
    # L = 258
    # lam = 0.5
    # m02s = np.linspace(-0.7, -0.6, 5 )


    # sweeps = 1000
    # cluster_method = WOLFF
    # thermalization = 20**2
    # record_rate = 10
    # cluster_rate = 5

    # quantities = ['magnetization', 'binder_cumulant', 'susceptibility']

    # data = {q : np.empty((2, len(m02s))) for q in quantities}

    # for i,m02 in enumerate(m02s):
        # if RANK==0: print('running '+str(i)+'...', end='\r')
        # start = time()
        # l = Phi4Lattice(dim=L, m02=m02, lam=lam)
        # rw = RandomWalk(l)

        # if RANK==0:
            # recorder = Recorder(gif=False)
        # else:
            # recorder = None
        # rw.run(sweeps, cluster_method=cluster_method, thermalization=thermalization, cluster_rate=cluster_rate, recorder=recorder, record_rate=record_rate)
        # if RANK==0:
            # print(f'{i} completed in {time() - start} seconds')
            # for q, arr in data.items():
                # arr[0,i] = recorder.mean(q)
                # arr[1,i] = recorder.error(q)

            # if recorder.gp is not None:
                # recorder.save_gif("converge.gif")

    # if RANK==0:
        # fig, axes = plt.subplots(3,1,figsize=(16,10))
        # for ax, (q, arr) in zip(axes, data.items()):
            # ax.errorbar(m02s, arr[0], arr[1], None, '.')
            # ax.set_ylabel(Recorder.quantities[q].label)
        # axes[-1].set_xlabel(r"$m_0^2$")
        # axes[0].set_title(f"Phase transition vs. $m_0^2$: L={L}, $\lambda={lam}$, 100 sweeps with Wolff every 5")

        # plt.show()


    # Test gradient flow
    L = 258
    lam = 0.5
    m02 = -0.7

    l = Phi4Lattice(dim=L, m02=m02, lam=lam)

    # sweeps = 1000
    # cluster_method = WOLFF
    # thermalization = 20**2
    # record_rate = 10
    # cluster_rate = 5

    quantities = ['magnetization', 'binder_cumulant', 'susceptibility', 'action']

    recorder = Recorder()
    taus = np.linspace(0.,0.005,100)

    for tau in taus:
        rho = deepcopy(l)
        rho.flow_evolve(tau)
        recorder.record(rho)

    fig, axes = plt.subplots(len(quantities), 1, figsize=(16, 3*len(quantities)))
    for ax, q in zip(axes, quantities):
        ax.plot(taus, recorder.values[q], '.')
        ax.set_ylabel(Recorder.quantities[q].label)
    axes[-1].set_xlabel(r"$\tau$")
    axes[0].set_title(f"Quantities in flow time; $L={L}$, $\\lambda={lam}$, $m_0^2={m02}$")

    plt.show()


    # just testing the gradient flow, so stop here

    # set the site assignments (factor of 2 for checkerboard)

if __name__=="__main__":
    main()
