
import numpy as np
import pickle
from matplotlib import pyplot as plt
from datetime import timedelta
import csv
from mcmc import Recorder
from gvar import gvar

L = 64
lam = 0.5
m02s = [-0.80,-0.76, -0.72, -0.68, -0.64]
o_filenames = [f"cpp_code/outputs/data_{m02}.csv" for m02 in m02s]

recorders = pickle.load(open('data/binder_cumulant.pickle', 'rb'))
quantities = ['magnetization','susceptibility','binder_cumulant']

print(len(recorders))
for r in recorders:
    r.finalize_values()

tds = []

plot_data_py = {
        'magnetizations' : [],
        'susceptibilities': [],
        'binder_cumulants': [],
}

def jackknife(f, phi, tau=0.5):
    cumsum = 0
    average = f(phi)
    for i in range(phi.size):
        cumsum += ( f(np.delete(phi,i)) - average )**2
    return np.sqrt(2 * tau * cumsum)

def magnetization(phi):
    return np.mean(np.abs(phi))

def binder_cumulant(phi):
    return 1 - np.mean(phi**4) / (3 * np.mean(phi**2)**2)

def susceptibility(phi):
    return np.mean(phi**2) - np.mean(np.abs(phi))**2

with open('data/log.out','r') as log_file:
    for line in log_file.readlines():
        if line[:8] == "Executed":
            time_str = line.split(" ")[-1]
            minutes,seconds = time_str.split(":")
            td = timedelta(minutes=int(minutes), seconds=int(seconds))
            tds.append(td.total_seconds()//60)

def bimodality(phi, dphi=0.1, **kwargs):
    bin_i = ((phi+0.5*dphi)//dphi).astype(int)
    reshuffled = np.where(bin_i<0, np.abs(bin_i), bin_i*2)
    bincounts = np.bincount(reshuffled)
    return 1 - bincounts[0] / np.amax(bincounts)

for f in [magnetization, susceptibility, binder_cumulant]:
    quantity = f.__name__

    means = []
    stds = []
    for r in recorders:
        phi = r.values['phi']
        means.append(f(phi))
        stds.append(jackknife(f, phi))

    # means = np.array([r.derived_values[quantity] for r in recorders])
    # stds  = np.array([r.derived_errors[quantity] for r in recorders])


    plot_data_py[quantity] = [gvar(r.derived_values[quantity], r.derived_errors[quantity]) for r in recorders]



#####################################################################



plot_data = {
        'magnetization' : [],
        'susceptibility': [],
        'binder_cumulant': [],
        'action' : []
}
labels = [r'$\langle|\bar\phi|\rangle$', r'$\chi$', r'$U$', r'$S$']

for m02, o_filenames in zip(m02s, o_filenames):
    actions = []
    phibars = []
    with open(o_filenames, 'r') as datafile:
        csv_data = csv.reader(datafile)
        for lines in csv_data:
            actions.append(float(lines[0]))
            phibars.append(float(lines[1]))
    actions = np.array(actions)
    phibars = np.array(phibars)

    for f in [magnetization, susceptibility, binder_cumulant]:
        plot_data[f.__name__].append( gvar(f(phibars), jackknife(f, phibars)) )

    plot_data['action'].append( gvar(0.,0.) )

def diff(x, y):
    # return 2 * (x - y) / (x + y)
    return (x - y)


fig, axes = plt.subplots(3,1,figsize=(16,10))

axes[0].set_title("Absolute difference of quantities between Python and C++ routines, 1000 measurements with 1000 sweep thermalization")
axes[0].set_ylabel(r"$\langle|\bar\phi|\rangle$")
axes[1].set_ylabel(r"$\chi$")
axes[-1].set_xlabel("$m_0^2$")

diff_plot_data = {}

for q in quantities:
    if q=="magnetization":
        print(plot_data[q][3], plot_data_py[q][3])
    diff_plot_data[q] = [diff(x, y) for x,y in zip(plot_data[q], plot_data_py[q]) ]

labels = [r'$\langle|\bar\phi|\rangle$', r'$\chi$', r'$U$', r'$S$']

for ax, d, l in zip(axes, diff_plot_data.values(), labels):
    ax.axhline(0.)
    ax.errorbar(m02s, [x.mean for x in d], yerr=[x.sdev for x in d], fmt='k.', capsize=3.)
    ax.set_ylabel(l)

plt.show()

