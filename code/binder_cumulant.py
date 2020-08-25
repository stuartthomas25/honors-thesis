from parallel_utils import RANK
from mcmc.keys import WOLFF, SWENDSEN_WANG
from mcmc import Phi4Lattice, RandomWalk, Recorder
import numpy as np
from time import time
import sys
if 'view' in sys.argv:
    from matplotlib import pyplot as plt
import pickle

SEED = True
Ls = [16, 32, 64]
lam = 0.5
# m02s = [-3, -2.5, -2, -1.5, -1]
m02s = [-0.80, -0.76, -0.72, -0.70, -0.68, -0.64]


# sweeps = 10**5
cluster_method = WOLFF
thermalization = 10**4
record_rate = 100
cluster_rate = 5

measurements = 10**4

sweeps = thermalization + record_rate * measurements



def main():

    recorders = []
    for m02 in m02s:
        if RANK==0: print(f"m02 = {m02:.2f}")
        for L in Ls:

            if RANK==0:
                recorder = Recorder(thermalization=thermalization, rate=record_rate, gif=False)
            else:
                recorder = None

            l = Phi4Lattice(dim=L, m02=m02, lam=lam)
            rw = RandomWalk(l)


            if RANK==0:
                start = time()

            rw.run(sweeps, cluster_method=cluster_method,  cluster_rate=cluster_rate, recorder=recorder, progress=True)

            if RANK==0:
                dt = time() - start
                print(f"Executed in {int(dt//60)}:{str(int(dt%60)).zfill(2)}")
                if recorder.gp is not None: recorder.save_gif("plots/BC.gif")

            recorders.append(recorder)

    if RANK==0:
        pickle.dump(recorders, open('data/binder_cumulant.pickle', 'wb'))
        print("Saved pickle")




def view():
    recorders = pickle.load(open('data/binder_cumulant.pickle', 'rb'))

    quantities = ['magnetization','susceptibility','binder_cumulant']

    # for i,r in enumerate(recorders):
        # n = r.record_count
        # assert n==101
        # phi2_avg = sum(r.values['phi2'])/n
        # phi4_avg = sum(r.values['phi4'])/n
        # plt.plot(r.values['phi4'], label='$\phi^4$')
        # plt.plot(r.values['phi2'], label='$\phi^2$')
        # plt.legend()
        # plt.show()
        # print(phi2_avg**2, phi4_avg)
        # print(i, 1 - phi4_avg/(3*phi2_avg**2))
        # print()

    # fig, axes = plt.subplots(4,1, figsize=(16,10))
    # for ax, quantity in zip(axes, Recorder.primary_observables):
        # for i,L in enumerate(Ls):
            # recorder = recorders[i+15]
            # values = recorder.values[quantity]
            # ax.plot(values, label=f"$N={L}$")
            # # ax.errorbar(m02s, derived_values, yerr=stds, label=f"$N={L}$")
            # ax.legend()

            # ax.set_ylabel(quantity)
    # plt.show()

    fig, axes = plt.subplots(3,1, figsize=(16,10))
    for ax, quantity in zip(axes, quantities):
        for i,L in enumerate(Ls):
            some_recorders = recorders[i::len(Ls)]
            derived_values = np.array([r.derived_value(quantity) for r in some_recorders])
            # stds  = np.array([r.errors[quantity] for r in some_recorders])
            stds  = np.array([0. for r in some_recorders]) # EDIT

            ax.errorbar(m02s, derived_values, yerr=stds, label=f"$N={L}$")
            ax.legend()
            if quantity=="binder_cumulant":
                ax.axhline(2/3, c='r')
            if quantity=="abs_magnetization":
                ax.set_ylim( (0.0, np.sqrt(-m02s[-1]/lam)) )
            # print("phi2",[r.values['phi2'] for r in some_recorders]) # EDIT

            ax.set_ylabel(Recorder.derived_observables[quantity].label)

    plt.show()

if __name__=="__main__":

    if 'view' in sys.argv:
        view()
    else:
        main()
