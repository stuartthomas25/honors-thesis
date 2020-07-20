from parallel_utils import RANK
from keys import WOLFF, SWENDSEN_WANG
from lattice import Phi4Lattice
from random_walk import RandomWalk
from matplotlib import pyplot as plt
from recorder import Recorder
import numpy as np
from time import time

def main():
    L = 258
    lam = 0.5
    m02s = np.linspace(-0.7, -0.6, 5 )


    sweeps = 1000
    cluster_method = WOLFF
    thermalization = 20**2
    record_rate = 10
    cluster_rate = 5

    quantities = ['magnetization', 'binder_cumulant', 'susceptibility']

    data = {q : np.empty((2, len(m02s))) for q in quantities}

    for i,m02 in enumerate(m02s):
        if RANK==0: print('running '+str(i)+'...', end='\r')
        start = time()
        l = Phi4Lattice(dim=L, m02=m02, lam=lam)
        rw = RandomWalk(l)

        if RANK==0:
            recorder = Recorder(gif=False)
        else:
            recorder = None

        rw.run(sweeps, cluster_method=cluster_method, thermalization=thermalization, cluster_rate=cluster_rate, recorder=recorder, record_rate=record_rate)
        if RANK==0:
            print(f'{i} completed in {time() - start} seconds')
            for q, arr in data.items():
                arr[0,i] = recorder.mean(q)
                arr[1,i] = recorder.error(q)

            if recorder.gp is not None:
                recorder.save_gif("converge.gif")

    if RANK==0:
        fig, axes = plt.subplots(3,1,figsize=(16,10))
        for ax, (q, arr) in zip(axes, data.items()):
            ax.errorbar(m02s, arr[0], arr[1], None, '.')
            ax.set_ylabel(Recorder.quantities[q].label)
        axes[-1].set_xlabel(r"$m_0^2$")
        axes[0].set_title(f"Phase transition vs. $m_0^2$: L={L}, $\lambda={lam}$, 100 sweeps with Wolff every 5")

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

    main()
