from parallel_utils import RANK
from keys import WOLFF, SWENDSEN_WANG
from lattice import Phi4Lattice
from random_walk import RandomWalk
import numpy as np
from recorder import Recorder
from time import time
import sys
if 'view' in sys.argv:
    from matplotlib import pyplot as plt
import pickle
from scipy.io import savemat
from UWerr import UWerr

L = 16
lam = 0.5
m02s = [-0.68]

# sweeps = 10**5
cluster_method = WOLFF
thermalization = 10**3
record_rate = 100
cluster_rate = 5

measurements = 10**3

sweeps = thermalization + record_rate * measurements
print(sweeps)
# print(sweeps, "sweeps")

quantities = ['magnetization', 'binder_cumulant', 'susceptibility']

def main():

    if RANK==0:
        recorder = Recorder(thermalization=thermalization, rate=record_rate, gif=False)
        data = np.empty((len(m02s), recorder.num_measurements(sweeps)))
    else:
        recorder = None


    for i, m02 in enumerate(m02s):

        l = Phi4Lattice(dim=L, m02=m02, lam=lam)
        rw = RandomWalk(l)
        start = time()
        rw.run(sweeps, cluster_method=cluster_method,  cluster_rate=cluster_rate, recorder=recorder, progress=True)

        if RANK==0:
            dt = time() - start
            print(f"Executed in {int(dt//60)}:{str(int(dt%60)).zfill(2)}")
            data[i,:] = recorder.values['magnetization']
            recorder.clear()

    if RANK==0:
        np.save('data/histogram', data)
        if recorder.gp is not None:
            recorder.save_gif('plots/lattice.gif')



def view():
    if RANK==0:
        data = np.load('data/histogram.npy')
        i = 0
        #UWerr(data[i]**2, name=f"$\\langle \\phi^2 \\rangle$: $L={L}$, $\\lambda={lam}$, $m_0^2={m02s[i]}$, {sweeps} sweeps, {thermalization} sweep thermalization, recording every {record_rate} sweeps", plot=True, whole_plot=True)
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
