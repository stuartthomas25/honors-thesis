from parallel_utils import RANK
from mcmc.keys import WOLFF, SWENDSEN_WANG
from mcmc import Phi4Lattice, RandomWalk, Recorder
from random import seed
import numpy as np
from time import time
import sys
if 'view' in sys.argv:
    from matplotlib import pyplot as plt
import pickle
from scipy.io import savemat
from UWerr import UWerr
import pickle

SEED = True
L = 64
lam = 0.5
m02s = [-0.72]

# sweeps = 10**5
cluster_method = None
thermalization = 10**4
record_rate = 1
cluster_rate = 5

measurements = 10**4

sweeps = thermalization + record_rate * measurements
# print(sweeps, "sweeps")

quantities = ['data', 'magnetization', 'susceptibility', 'binder_cumulant', 'action']

def main():

    if RANK==0:
        recorders = [Recorder(quantities, thermalization=thermalization, rate=record_rate, gif=False) for _ in m02s]
    else:
        recorders = [None for _ in m02s]


    for i, (recorder, m02) in enumerate(zip(recorders, m02s)):

        l = Phi4Lattice(dim=L, m02=m02, lam=lam)
        # l.data[:,:] = np.sqrt(-m02/lam)
        # l.action = l.calculate_action()

        rw = RandomWalk(l)
        start = time()
        if RANK==0: recorder.record(l)
        rw.run(sweeps, cluster_method=cluster_method,  cluster_rate=cluster_rate, recorder=recorder, progress=True)

        if RANK==0:
            dt = time() - start
            print(f"Executed in {int(dt//60)}:{str(int(dt%60)).zfill(2)}")

    if RANK==0:
        rw.check_action()

        recorders[i].record(l)
        pickle.dump(recorders, open('data/converge.pickle', 'wb'))
        print("Saved pickle")
        if recorders[i].gp is not None:
            recorders[i].save_gif('plots/lattice.gif')
            print("Saved gif")



def view():
    if RANK==0:
        recorders = pickle.load(open('data/converge.pickle', 'rb'))
        i = 0
        recorder=recorders[i]
        data = np.array(recorder.values['data'])
        phi_sq = np.mean(data**2, axis=(1,2))

        UWerr(phi_sq, name=f"$\\langle \\phi^2 \\rangle$: $L={L}$, $\\lambda={lam}$, $m_0^2={m02s[i]}$, {sweeps} sweeps, {thermalization} sweep thermalization, recording every {record_rate} sweeps", plot=True, whole_plot=True)

        fig, axes = plt.subplots(4, 1, figsize=(16,16))
        sweep_x = np.arange(recorder.record_count)

        phi_broken =np.sqrt(-m02s[i]/lam)

        axes[0].axhline(phi_broken, ls=':', c='k')
        axes[0].axhline(0, ls=':', c='k')
        axes[0].plot(sweep_x, np.abs(recorder.values['magnetization']))
        axes[0].set_ylabel(r'$\langle \phi \rangle$')


        axes[1].set_ylabel("$\chi$")
        axes[1].plot(sweep_x, recorder.values['susceptibility'])
        axes[2].set_ylabel("$U$")
        axes[2].plot(sweep_x, recorder.values['binder_cumulant'])
        axes[3].set_ylabel("$S$")
        axes[3].plot(sweep_x, recorder.values['action'])
        axes[3].axhline(0., ls=':', c='r')
        axes[3].axhline(-0.25 * m02s[i]**2/lam, ls=':', c='k')

        plt.show()
        plt.hist(data[-1].flatten(), bins=100)
        plt.axvline(-np.sqrt(-m02s[i]/lam),c='k')
        plt.axvline( np.sqrt(-m02s[i]/lam),c='k')
        plt.show()
        quit()

        nseries = data.shape[0]
        fig, axes = plt.subplots(nseries, 1, sharex=True, figsize=(16, 8*nseries), squeeze=False)
        axes = [ax for row in axes for ax in row] # flatten
        for m02, series, ax in zip(m02s, data, axes):
            ax.hist(series, 100, label=f'$m_0^2={m02}$', alpha=0.7 )
            ax.legend()
            # ax.set_xlim((-0.75,0.75))
        axes[0].set_title(f"$\\langle \\phi \\rangle$ histograms: $L={L}$, $\\lambda={lam}$, {sweeps} sweeps, {thermalization} sweep thermalization, recording every {record_rate} sweeps")
        axes[-1].set_xlabel(r"$\langle \phi \rangle$")
        plt.show()


# # Test gradient flow
    # recorder = Recorder()

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

    # # set the site assignments (factor of 2 for checkerboard)


    # cluster_rate = 1
    # record_rate = 1
    # thermalization = 1

    # recorder = Recorder()
    # recorder.save(l)

    # rw = RandomWalk(l)

    # start = time()

    # rw.run(sweeps, cluster_method, cluster_rate, record_rate, thermalization, recorder)

    # if RANK==0:
        # exec_time = time() - start
        # print(f"{exec_time}")

        # recorder.save_gif("plots/lattice_visualization.gif")

        # recorder.plot(title=f"Monte Carlo Simulation of $\phi^4$ Model using Metropolis and Wolff Algorithms, $L={L}$, $\\lambda={lam}$, $\\mu_0^2={m02}$, $t={exec_time:.1f}s$", fname="parallel.png", show=('-S' in dopts))

if __name__=="__main__":
    if 'view' in sys.argv:
        view()
    else:
        main()
