from mpi4py import MPI
from matplotlib import pyplot as plt
import imageio

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()

class Recorder(object):
    def __init__(self, gif = False):
        if RANK==0:
            if gif:
                self.gp = GifProducer()
            else:
                self.gp = None
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
            self.magnetizations.append(l.magnetization())
            self.energies.append(l.action/l.size)
            if self.gp is not None:
                self.gp.save_lat(l)
            self.record_count += 1

    def save_gif(self, fname, fps=3):
        if RANK==0:
            self.gp.save(fname, fps)

    def plot(self, title="", fname=None, show=True, xlabel="Sweep"):
        if RANK==0:
            if fname is None:
                fname = "Figure.png"
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
            ax4.set_ylim((-0.1,1.1))
            ax4.set_xlabel(xlabel)
            # ax3.set_ylim((-0.05,1.05))
            ax1.set_title(title)
            plt.savefig('plots/'+fname)
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
