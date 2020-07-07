from mpi4py import MPI
from matplotlib import pyplot as plt
import imageio

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()

class Recorder(object):
    def __init__(self):
        if RANK==0:
            self.gp = GifProducer()
            self.iter_index = []
            self.magnetizations = []
            self.energies = []
            self.susceptibilities = []
            self.binder_cums = []
            self.states = []
            self.record_count = 0
            self.thermalization = 50

    def save(self, l):
        if RANK==0:
            self.iter_index.append(self.record_count)
            self.susceptibilities.append(l.susceptibility())
            self.binder_cums.append(l.binder_cumulant())
            self.magnetizations.append(abs(l.magnetization()))
            l.calculate_action()
            self.energies.append(l.action/l.size)
            self.gp.save_lat(l)
            self.record_count += 1

    def save_gif(self, fname):
        if RANK==0:
            self.gp.save(fname)

    def plot(self, L, m, lam, exec_time, show):
        if RANK==0:
            x = self.iter_index
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(18,12), sharex = True)
            ax1.plot(x, self.magnetizations)
            ax2.plot(x, self.energies)
            ax3.plot(x, self.susceptibilities)
            ax4.plot(x, self.binder_cums)

            ax1.set_ylabel("Magnetization")
            ax2.set_ylabel("Total Action")
            ax3.set_ylabel("Susceptibility")
            ax4.set_ylabel("Binder Cumulant")
            ax4.set_xlabel("Sweep")
            # ax3.set_ylim((-0.05,1.05))
            ax1.set_title(f"Monte Carlo Simulation of $\phi^4$ Model using Metropolis and Wolff Algorithms, $L={L}$, $\\lambda={lam}$, $\\mu_0^2={m**2}$, $t={exec_time:.1f}s$")
            plt.savefig('plots/temp_parallel.png')
            if show:
                plt.show()

class GifProducer(object):
    def __init__(self):
        self.frames = []

    def save_lat(self, lat):
        im = lat.show(show=False)
        self.frames.append(im.make_image("AGG")[0])
        plt.close()

    def save(self, fn, fps=3):
        imageio.mimwrite(fn, self.frames, fps=fps)
