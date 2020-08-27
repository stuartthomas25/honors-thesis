from mpi4py import MPI
from matplotlib import pyplot as plt
import imageio
from functools import wraps
from math import sqrt
import numpy as np
import gvar

# COMM = MPI.COMM_WORLD
# RANK = COMM.Get_rank()


class Quantity(object):
    def __init__(self, f, label=None):
        self.f = f
        self.label = label

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)

class Recorder(object):

    primary_observables = {} # phi
    secondary_observables = {} # |phi|, phi^2, phi^4, ...
    derived_observables = {} # U, chi, mu, etc

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

    @classmethod
    def secondary_observable(cls, f):
        name = f.__name__
        cls.secondary_observables[name] = Quantity(f)
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

    @property
    def gvars(self):
        primary_keys = []
        secondary_keys = []

        all_values = np.empty((len(type(self).primary_observables) + len(type(self).secondary_observables), self.record_count))

        for i,(k,pv) in enumerate(self.values.items()):
            primary_keys.append(k)
            all_values[i,:] = pv
        for i,(k,q) in enumerate(type(self).secondary_observables.items()):
            offset = len(primary_keys)
            secondary_keys.append(k)
            for j in range(self.record_count):
                sv = q(**{key:val[j] for  key,val in self.values.items()})
                all_values[i+offset,j] = sv

        cov_mat = np.cov(all_values)
        means = np.mean(all_values, axis=1)


        gvars = gvar.gvar(means, cov_mat, verify=True)
#         print(*[gv.sdev for gv in gvars], sep='\t')
        # print(*[np.std(arr)/np.sqrt(self.record_count) for arr in all_values], sep='\t')
        # print(*[np.std(arr) for arr in all_values], sep='\t')
        # print()


        return {k:gv for k,gv in zip(primary_keys + secondary_keys, gvars)}


        # for k, q in type(self).secondary_observables.items():
            # secondary_values[k] = np.array([q(**{k:v[i] for k,v in self.values.items()}) for i in range(self.record_count)])
        # obsrv_matrix

            # errors = np.std(nparr) / sqrt(nparr.size)
        # return {k : gvar.gvar(means[k], errors[k]) for k in self.values}

    def derived_value(self, key):
        return type(self).derived_observables[key] (**self.means)

    def derived_gvar(self, key):
        return type(self).derived_observables[key] (**self.gvars)

@Recorder.primary_observable
def phi(lat):
    return np.sum(lat.data) / lat.size

@Recorder.secondary_observable
def abs_phi(phi):
    return abs(phi)

@Recorder.secondary_observable
def phi2(phi):
    return phi**2

@Recorder.secondary_observable
def phi4(phi):
    return phi**4




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

