from mpi4py import MPI
# from matplotlib import pyplot as plt
import imageio
from functools import wraps
from math import sqrt
import numpy as np
from math import sqrt

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


    __statevars__ = ['values', 'record_count', 'rate', 'thermalization']

    def __init__(self, rate=1, thermalization=0, gif = False):
        if gif:
            self.gp = GifProducer()
        else:
            self.gp = None

        self.thermalization = thermalization
        self.rate = rate

        self.values = {}

        self.means = {}
        self.derived_values = {}
        self.derived_errors = {}

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

    def jackknife_means(self, i):
        ''' i indicates measurement to leave out (for jackknife)'''
        means = {}
        for k, arr in self.values.items():
            means[k] = (np.sum(arr) - arr[i]) / (self.record_count-1)
        return means

    def finalize_values(self, tau=0.5):
        # make all lists numpy arrays
        for k in self.values:
            self.values[k] = np.array(self.values[k])

        # add secondary observables
        all_values = dict(self.values)
        for k,q in Recorder.secondary_observables.items():
            all_values[k] = np.array(q(**self.values))
        self.values = all_values

        # calculate primary and secondary means
        for k, arr in self.values.items():
            self.means[k] = np.sum(self.values[k]) / self.record_count

        for k,q in Recorder.derived_observables.items():
            self.derived_values[k] = q(**self.means)

        cumsums = {k:0. for k in Recorder.derived_observables}
        for i in range(self.record_count):
            means = self.jackknife_means(i)
            for k,q in Recorder.derived_observables.items():
                cumsums[k] += (q(**means) - self.derived_values[k]) **2

        for k in cumsums:
            self.derived_errors[k] = sqrt(2 * tau * cumsums[k])

#     def gvars(self):
        # primary_keys = []
        # secondary_keys = []

        # all_values = np.empty((len(type(self).primary_observables) + len(type(self).secondary_observables), self.record_count))

        # for i,(k,pv) in enumerate(self.values.items()):
            # primary_keys.append(k)
            # all_values[i,:] = pv
        # for i,(k,q) in enumerate(type(self).secondary_observables.items()):
            # offset = len(primary_keys)
            # secondary_keys.append(k)
            # for j in range(self.record_count):
                # sv = q(**{key:val[j] for  key,val in self.values.items()})
                # all_values[i+offset,j] = sv

        # cov_mat = np.cov(all_values)
        # means = np.mean(all_values, axis=1)


        # gvars = gvar.gvar(means, cov_mat, verify=True)
# #         print(*[gv.sdev for gv in gvars], sep='\t')
        # # print(*[np.std(arr)/np.sqrt(self.record_count) for arr in all_values], sep='\t')
        # # print(*[np.std(arr) for arr in all_values], sep='\t')
        # # print()


        # return {k:gv for k,gv in zip(primary_keys + secondary_keys, gvars)}


        # # for k, q in type(self).secondary_observables.items():
            # # secondary_values[k] = np.array([q(**{k:v[i] for k,v in self.values.items()}) for i in range(self.record_count)])
        # # obsrv_matrix

            # # errors = np.std(nparr) / sqrt(nparr.size)
        # # return {k : gvar.gvar(means[k], errors[k]) for k in self.values}

    def derived_value(self, key):
        return type(self).derived_observables[key] (**self.means())

    def derived_gvar(self, key):
        return type(self).derived_observables[key] (**self.gvars)

    def derived_error(self, key, tau=0.5):
        cumsum = 0.
        mu = self.derived_value(key)
        for i in range(self.record_count):
            cumsum += (type(self).derived_observables[key] (**self.means(i)) - mu) **2

        return sqrt(2 * tau * cumsum)

    def __getstate__(self):
        d = {}
        for v in Recorder.__statevars__:
            d[v] = self.__dict__[v]
        return d

    def __setstate__(self, d):
        self.__init__()
        for k,v in d.items():
            self.__dict__[k] = v

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

