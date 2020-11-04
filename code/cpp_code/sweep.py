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

actions = []
phibars = []
# m02s = [-0.80, -0.76, -0.72, -0.68, -0.64];
# m02s = np.linspace(-0.80,-0.64, 17)
# m02s = np.linspace(-1.80, 1.00, 29)
m02s = np.linspace(-0.8, -0.64, 9)
# m02s = [-0.80]
o_filenames = [f"outputs/data_{m02}.csv" for m02 in m02s]
# streams = [StringIO for _ in m02s]
processes = []
L = 32

refresh_rate = 0.1

# some functions

def statvar(arr, tau=0.5):
    arr = np.array(arr)
    mean = np.mean(arr)

    # use jackknife to calculate errs
    def C_i(arr, i=None):
        ''' i indicates measurement to leave out (for jackknife)'''
        return (np.sum(arr) - arr[i]) / (arr.size-1)

    cumsum = 0
    for i in range(arr.size):
       cumsum += (C_i(arr, i) - mean)**2
    stddev = np.sqrt(2 * tau * cumsum)
    return gvar(mean, stddev)


def magnetization(phi):
    return statvar(np.abs(phi))

def binder_cumulant(phi):
    return 1 - statvar(phi**4) / (3 * statvar(phi**2)**2)

def susceptibility(phi):
    return statvar(phi**2) - statvar(np.abs(phi))**2

if "run" in sys.argv:
    quiet = "-q" in sys.argv
    if not quiet:
        term = Terminal()
        if term.height+1<len(m02s):
            raise Exception("terminal too short")
    for i,(m02, o_filename) in enumerate(zip(m02s, o_filenames)):
        cmd = ["./bin/sweep", "-m", str(m02), "-L", str(L), "-o", o_filename]
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
    fig, axes = plt.subplots(4,1,figsize=(16,10))

    axes[0].set_ylabel(r"$\langle|\bar\phi|\rangle$")
    axes[1].set_ylabel(r"$\chi$")
    axes[2].set_ylabel(r"$U$")
    axes[3].set_ylabel("S")
    axes[-1].set_xlabel("$m_0^2$")

    plot_data = {
            'magnetizations' : [],
            'susceptibilities': [],
            'binder_cumulants': [],
            'actions' : []
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

        plot_data['magnetizations'].append( magnetization(phibars) )
        plot_data['susceptibilities'].append( susceptibility(phibars) )
        plot_data['binder_cumulants'].append( binder_cumulant(phibars) )
        plot_data['actions'].append( statvar(actions) )

    for ax, d, l in zip(axes, plot_data.values(), labels):
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
