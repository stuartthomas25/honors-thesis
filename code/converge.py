from parallel_utils import RANK
from keys import WOLFF, SWENDSEN_WANG
from lattice import Phi4Lattice
from random_walk import RandomWalk
from matplotlib import pyplot as plt
from recorder import Recorder
import numpy as np
from time import time

def main():
    L = 64
    lam = 0.5
    m02s = np.linspace(-1., -0.5, 30)


    sweeps = 1000
    cluster_method = WOLFF
    thermalization = 20**2
    record_rate = 100
    cluster_rate = 5


    measurements = (sweeps-thermalization) // record_rate - 1
    m = np.empty((len(m02s), measurements))
    bc = np.empty((len(m02s), measurements))
    s = np.empty((len(m02s), measurements))

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
            m[i,:] = recorder.magnetizations
            bc[i,:] = recorder.binder_cums
            s[i,:] = recorder.susceptibilities
            print(f'{i} completed in {time() - start} seconds')

            if recorder.gp is not None:
                recorder.save_gif("converge.gif")

    if RANK==0:
        fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(16,10))
        ax1.plot(m02s, m, '.')
        ax2.plot(m02s, bc, '.')
        ax3.plot(m02s, s, '.')
        ax1.set_ylabel(r"$\langle \phi \rangle$")
        ax2.set_ylabel(r"BC")
        ax3.set_ylabel(r"S")
        ax3.set_xlabel(r"$m_0^2$")
        ax1.set_title(f"Phase transition vs. $m_0^2$: L={L}, $\lambda={lam}$, 100 sweeps with Wolff every 5")

        plt.show()


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
