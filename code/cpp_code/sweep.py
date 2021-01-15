import numpy as np

from matplotlib import pyplot as plt
import csv
import os
import sys
from blessings import Terminal
from subprocess import Popen, PIPE, CalledProcessError
# from io import StringIO
import time
from gvar import gvar

cores = 1
actions = []
phibars = []
# m02s = [-0.80, -0.76, -0.72, -0.68, -0.64];
# m02s = np.linspace(-0.80,-0.64, 17)
# m02s = np.linspace(-1.80, 1.00, 29)
m02s = [-0.80,-0.76, -0.72, -0.68, -0.64]
# m02s = [-0.80]
o_filenames = [f"outputs/data_{m02}.csv" for m02 in m02s]
# streams = [StringIO for _ in m02s]
processes = []
L = 64

refresh_rate = 0.1

# some functions


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

if "run" in sys.argv:
    quiet = "-q" in sys.argv
    if not quiet:
        term = Terminal()
        if term.height+1<len(m02s):
            raise Exception("terminal too short")
    for i,(m02, o_filename) in enumerate(zip(m02s, o_filenames)):
        cmd = ["mpiexec","-n", str(cores), "bin/sweep", "-m", str(m02), "-L", str(L), "-o", o_filename]
        # cmd = ["./bin/sweep", "-m", str(m02), "-L", str(L), "-o", o_filename]
        processes.append( Popen(cmd, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True) )

    #Now watch the results
    if not quiet:
        for _ in processes:
            print("...")

    started = [False for _ in processes]
    finished = [False for _ in processes]
    while not all(finished):
        for i,p in enumerate(processes):

            line = p.stdout.readline()[:-1]
            if line != "":
                if not started[i]:
                    started[i] = True
                if not quiet:
                    terminal_line = term.height-len(processes)-1+i
                    with term.location(0,terminal_line):
                        print( f" {i}> " + line)
            else:
                if started[i]:
                    finished[i] = True





elif "view" in sys.argv:
    fig, axes = plt.subplots(3,1,figsize=(16,10), sharex=True)

    axes[0].set_ylabel(r"$\langle|\bar\phi|\rangle$")
    axes[1].set_ylabel(r"$\chi$")
    axes[2].set_ylabel(r"$U$")
    # axes[3].set_ylabel("S")
    axes[-1].set_xlabel("$m_0^2$")

    plot_data = {
            'magnetization' : [],
            'susceptibility': [],
            'binder_cumulant': [],
            # 'action' : []
    }
    labels = [r'$\langle|\bar\phi|\rangle$', r'$\chi$', r'$U$']

    for m02, o_filename in zip(m02s, o_filenames):
        actions = []
        phibars = []
        with open(o_filename, 'r') as datafile:
            csv_data = csv.reader(datafile)
            for lines in csv_data:
                actions.append(float(lines[0]))
                phibars.append(float(lines[1]))
        actions = np.array(actions)
        phibars = np.array(phibars)

        for f in [magnetization, susceptibility, binder_cumulant]:
            plot_data[f.__name__].append( gvar(f(phibars), jackknife(f, phibars)) )

        # plot_data['action'].append( gvar(0.,0.) )

    for ax, d, l in zip(axes, plot_data.values(), labels):
        if ax==axes[1]:
            print(d[2].mean, d[2].sdev)
        ax.errorbar(m02s, [x.mean for x in d], yerr=[x.sdev for x in d], fmt='k.', capsize=3.)
        ax.set_ylabel(l)

    axes[0].set_ylim((0., 1.1))
    axes[2].axhline(2/3,ls=':')
    axes[2].set_ylim((-0.1, None))
    plt.show()

elif "evolution" in sys.argv:
    fig, axes = plt.subplots(2,1,figsize=(8,10))

    axes[0].set_ylabel(r"$\bar\phi$")
    axes[1].set_ylabel("S")
    axes[-1].set_xlabel("sweeps")
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

        axes[0].plot(np.abs(phibars), '-', label=f"$m_0^2$={m02}")
        axes[1].plot(actions, label=f"$m_0^2$={m02}")

    plt.legend()
    plt.show()

else:
    print("No command given, please use 'run' or 'view' arguments")
