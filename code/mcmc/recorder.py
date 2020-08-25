from mpi4py import MPI
from matplotlib import pyplot as plt
import imageio
from functools import wraps
from math import sqrt
import numpy as np

# COMM = MPI.COMM_WORLD
# RANK = COMM.Get_rank()


class Quantity(object):
    def __init__(self, f, label=None):
        self.f = f
        self.label = label

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)

class Recorder(object):

    primary_observables = {}
    derived_observables = {}

    def __init__(self, rate=1, thermalization=0, gif = False):
        if gif:
            self.gp = GifProducer()
        else:
            self.gp = None

        self.thermalization = thermalization
        self.rate = rate

        self.values = {}

        for obsrv in type(self).primary_observables:
            self.values[obsrv] = []

        self.record_count = 0

        # print("Recorder",Recorder.quantities)

    @classmethod
    def derived_observable(cls, label):
        def decorator(f):
            name = f.__name__
            cls.derived_observables[name] = Quantity(f, label)
            return f
        return decorator

    @classmethod
    def primary_observable(cls, f):
        name = f.__name__
        cls.primary_observables[name] = Quantity(f)
        return f


    def num_measurements(self, sweeps):
        return (sweeps - self.thermalization) // self.rate - 1

    def clear(self):
        for v in self.values:
            self.values[v] = []

    def record(self, lat):
        for val, arr in self.values.items():
            arr.append(type(self).primary_observables[val] (lat))
        if self.gp is not None:
            self.gp.save_lat(lat)
        self.record_count += 1

    def save_gif(self, fname, fps=3):
        self.gp.save(fname, fps)

    @property
    def means(self):
        return {k : sum(self.values[k]) / self.record_count for k in self.values}

    @property
    def errors(self):
        stderrs = {}
        for key in self.values:
            nparr = np.array(self.values[key])
            stddev = np.std(nparr)
            stderrs[key] = stddev / sqrt(nparr.size)

        return stderrs

    def derived_value(self, key):
        return type(self).derived_observables[key] (**self.means)

@Recorder.primary_observable
def phi(lat):
    return np.sum(lat.data) / lat.size

@Recorder.primary_observable
def abs_phi(lat):
    return abs(np.sum(lat.data)) / lat.size

@Recorder.primary_observable
def phi2(lat):
    # return np.sum(lat.data**2) / lat.size
    return (np.sum(lat.data) / lat.size) **2

@Recorder.primary_observable
def phi4(lat):
    return (np.sum(lat.data) / lat.size) **4




@Recorder.derived_observable(r"$\langle|\bar \phi|\rangle$")
def magnetization(abs_phi, **kwargs):
    return abs_phi

@Recorder.derived_observable(r"$U$")
def binder_cumulant(phi4, phi2, **kwargs):
    return 1 - phi4 / (3 * phi2**2)

@Recorder.derived_observable(r"$\chi$")
def susceptibility(phi2, abs_phi, **kwargs):
    return phi2 - abs_phi**2


class GifProducer(object):
    def __init__(self):
        self.frames = []

    def save_lat(self, lat):
        im = lat.show(show=False)
        self.frames.append(im.make_image("AGG")[0])
        plt.close()

    def save(self, fn, fps=3):
        imageio.mimwrite(fn, self.frames, fps=fps)

