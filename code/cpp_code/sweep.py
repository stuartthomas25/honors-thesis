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
betas = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
# m02s = [-0.80]
o_filenames = [f"outputs/data_{beta}.csv" for beta in betas]
# streams = [StringIO for _ in m02s]
processes = []
L = 128

refresh_rate = 0.1

# some functions


def jackknife(f, phi, S, beta, tau=0.5):
    cumsum = 0
    average = f(phi,S,beta)
    for i in range(phi.size):
        cumsum += ( f(np.delete(phi,i), np.delete(S,i),beta) - average )**2
    return np.sqrt(2 * tau * cumsum)

def magnetization(phi):
    return np.mean(np.abs(phi))

def binder_cumulant(phi):
    return 1 - np.mean(phi**4) / (3 * np.mean(phi**2)**2)

def susceptibility(phi):
    return np.mean(phi**2) - np.mean(np.abs(phi))**2

def action(phi,S,beta):
    return np.mean(S)/(beta*L**2)

def internal_energy(phi,S,beta):
    return np.mean(S)/(0.5*beta*L**2)

if "run" in sys.argv:
    quiet = "-q" in sys.argv
    if not quiet:
        term = Terminal()
        if term.height+1<len(betas):
            raise Exception("terminal too short")
    for i,(beta, o_filename) in enumerate(zip(betas, o_filenames)):
        cmd = ["mpiexec","-n", str(cores), "bin/sweep", "-b", str(beta), "-L", str(L), "-o", o_filename]
        if quiet: cmd.append('-q')
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
            if "Working" in line:
                if not started[i]:
                    started[i] = True
                if not quiet:
                    terminal_line = term.height-len(processes)-1+i
                    with term.location(0,terminal_line):
                        print( f" {i}> " + line)
            else:
                line = p.stdout.readline()[:-1]
                if started[i]:
                    finished[i] = True
                if not quiet:
                    terminal_line = term.height-len(processes)-1+i
                    with term.location(0,terminal_line):
                        print( f" {i}| " + line)



elif "view" in sys.argv:
    fig, axes = plt.subplots(1,1,figsize=(16,10), sharex=True, squeeze=False)
    axes = [ax for row in axes for ax in row]

    axes[0].set_ylabel(r"$E$")
    axes[0].set_xlabel(r"$\beta$")

    plot_data = {
            'internal_energy' : []
    }
    labels = ["$S$"]

    for beta, o_filename in zip(betas, o_filenames):
        actions = []
        phibars = []
        with open(o_filename, 'r') as datafile:
            csv_data = csv.reader(datafile)
            for lines in csv_data:
                actions.append(float(lines[0]))
                phibars.append(float(lines[1]))
        actions = np.array(actions)
        phibars = np.array(phibars)


        for f in [internal_energy]:
            plot_data[f.__name__].append( gvar(f(phibars, actions, beta), jackknife(f, phibars, actions, beta)) )

        # plot_data['action'].append( gvar(0.,0.) )

    for ax, d, l in zip(axes, plot_data.values(), labels):
        ax.errorbar(betas, [x.mean for x in d], yerr=[x.sdev for x in d], fmt='k.', capsize=3.)
        ax.set_ylabel(l)



    # Theoretical values
    beta_weak = np.arange(0.001,1.7, 0.1)
    beta_strong = np.arange(1.2, 3.4, 0.1)

    y = np.cosh(beta_weak)/np.sinh(beta_weak) - 1/beta_weak
    Eweak = 4 - 4*y - 8 * y**3 - 48/5 * y**5

    Estrong = 2/beta_strong + 1/(4*beta_strong**2) + 0.156/beta_strong**3


    axes[0].plot(beta_weak, Eweak)
    axes[0].plot(beta_strong, Estrong)

    # axes[0].set_ylim((0., 4.))
    plt.show()

elif "evolution" in sys.argv:
    fig, axes = plt.subplots(2,1,figsize=(8,10))

    axes[0].set_ylabel(r"$\bar\phi$")
    axes[1].set_ylabel("S")
    axes[-1].set_xlabel("sweeps")
    for beta, o_filenames in zip(betas, o_filenames):
        actions = []
        phibars = []
        with open(o_filenames, 'r') as datafile:
            csv_data = csv.reader(datafile)
            for lines in csv_data:
                actions.append(float(lines[0]))
                phibars.append(float(lines[1]))

        actions = np.array(actions)
        phibars = np.array(phibars)

        axes[0].plot(np.abs(phibars), '-', label=f"$\\beta^2$={beta}")
        axes[1].plot(actions, label=f"$\\beta$={beta}")

    plt.legend()
    plt.show()

else:
    print("No command given, please use 'run' or 'view' arguments")
