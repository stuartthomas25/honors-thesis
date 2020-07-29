from parallel_utils import RANK
from keys import WOLFF, SWENDSEN_WANG
from lattice import Phi4Lattice
from random_walk import RandomWalk
from matplotlib import pyplot as plt
from recorder import Recorder
import numpy as np
from time import time
import sys
import pickle
from scipy.io import savemat
from UWerr import UWerr

Ls = [16, 32, 64, 128]
lam = 0.5
m02 = -0.72

cluster_method = WOLFF
thermalization = 10**4
cluster_rate = 5

def main():

    if RANK==0:
        datas = {}

    for L in Ls:

        l = Phi4Lattice(dim=L, m02=m02, lam=lam)

#         phis = np.linspace(-2,2,200)
        # Vs = np.empty_like(phis)

        # for i,phi in enumerate(phis):
           # Vs[i] = l.static_lagrangian(phi, [phi, phi])

        # plt.plot(phis, Vs)
        # plt.axvline(np.sqrt(-m02/lam))
        # plt.axvline(-np.sqrt(-m02/lam))

        # plt.show()


        # quit()

        if RANK==0:
            recorder = Recorder(thermalization=0, rate=10, gif=False)
        else:
            recorder = None



        rw = RandomWalk(l)
        start = time()
        rw.run(thermalization, cluster_method=cluster_method,  cluster_rate=cluster_rate, recorder=recorder, progress=True)

        if RANK==0:
            print(f"Executed in {(time() - start)} sec")
            datas[L] = l.data.flatten().tolist()
            # recorder.save_gif('plots/lattice.gif')

    if RANK==0:
        pickle.dump(datas, open("data/histogram.pickle", 'wb'))

def view():
    if RANK==0:
        datas = pickle.load(open('data/histogram.pickle', 'rb'))

        fig, ax = plt.subplots(1, 1, figsize=(16, 8))


        for L, data in reversed(datas.items()):
            ax.hist(data, int(L), label=f'$L={L}$')

        ax.axvline( np.sqrt(-m02/lam),c='k')
        ax.axvline(-np.sqrt(-m02/lam),c='k')
        ax.legend()
        ax.set_xlim((-3,3))

        ax.set_title(f"$\\langle \\phi \\rangle$ histograms: $L={L}$, $\\lambda={lam}$, $m_0^2={m02}$, {thermalization} sweep thermalization")
        ax.set_xlabel(r"$\langle \phi \rangle$")
        plt.show()

if __name__=="__main__":
    if 'view' in sys.argv:
        view()
    else:
        main()
