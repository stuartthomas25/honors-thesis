from mpi4py import MPI
from matplotlib import pyplot as plt
import imageio
from functools import wraps
from math import sqrt
import numpy as np

# COMM = MPI.COMM_WORLD
# RANK = COMM.Get_rank()


class Quantity(object):
    def __init__(self, f, label):
        self.f = f
        self.label = label

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)


class Recorder(object):

    quantities = {}

    def __init__(self, gif = False):
        if gif:
            self.gp = GifProducer()
        else:
            self.gp = None

        self.xvals = []
        self.values = {}
        for quant in type(self).quantities:
            self.values[quant] = []

        self.record_count = 0

        # print("Recorder",Recorder.quantities)

    @classmethod
    def quantity(cls, label):
        def decorator(f):
            name = f.__name__
            cls.quantities[name] = Quantity(f, label)
            return f
        return decorator


    def record(self, lat):
        self.xvals.append(self.record_count)
        for val, arr in self.values.items():
            arr.append(type(self).quantities[val] (lat))
        if self.gp is not None:
            self.gp.save_lat(lat)
        self.record_count += 1

    def save_gif(self, fname, fps=3):
        self.gp.save(fname, fps)

    # def plot(self, title="", fname=None, show=True, xlabel="Sweep"):
        # if fname is None:
            # fname = "Figure.png"
        # x = self.iter_index
        # fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(18,12), sharex = True)
        # ax1.plot(x, self.magnetizations)
        # ax2.plot(x, self.energies)
        # ax3.plot(x, self.susceptibilities)
        # ax4.plot(x, self.binder_cums)

        # ax1.set_ylabel("Magnetization")
        # ax2.set_ylabel("Total Action")
        # ax3.set_ylabel("Susceptibility")
        # ax4.set_ylabel("Binder Cumulant")
        # ax4.set_ylim((-0.1,1.1))
        # ax4.set_xlabel(xlabel)
        # # ax3.set_ylim((-0.05,1.05))
        # ax1.set_title(title)
        # plt.savefig('plots/'+fname)
        # if show:
            # plt.show()

    def mean(self, key):
        arr = self.values[key]
        return sum(arr) / len(arr)

    def error(self, key):
        if key not in self.values:
            raise Exception("Quantity key not found")
        nparr = np.array(self.values[key])
        stddev = np.std(nparr)
        stderr = stddev / sqrt(nparr.size)

        return stderr

@Recorder.quantity(r"$|\langle \phi \rangle|$")
def magnetization(lat):
    return abs(np.sum(lat.data) / lat.size)

@Recorder.quantity(r"$U$")
def binder_cumulant(lat):
    phi_sq = np.sum(lat.data**2) / lat.size
    phi_qu = np.sum(lat.data**4) / lat.size
    return 1 - phi_qu / (3 * phi_sq**2)

@Recorder.quantity(r"$\chi$")
def susceptibility(lat):
    m =  np.sum(lat.data) / lat.size
    m2 = np.sum(lat.data**2) / lat.size
    return (m2 - m**2)

@Recorder.quantity(r"$S$")
def action(lat):
    return lat.action / lat.size

class GifProducer(object):
    def __init__(self):
        self.frames = []

    def save_lat(self, lat):
        im = lat.show(show=False)
        self.frames.append(im.make_image("AGG")[0])
        plt.close()

    def save(self, fn, fps=3):
        imageio.mimwrite(fn, self.frames, fps=fps)
