#!/opt/miniconda3/envs/thesis/bin python

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
from inspect import signature
import click
from UWerr import UWerr
import yaml

delay = 0.4
cores = 1
actions = []
phibars = []
# betas = [0.6, 0.8, 1.0,1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
orig_betas = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
betas = [1.0, 0.5, 2.5, 1.55, 1.65, 1.75]
betas = betas+orig_betas
# betas = [1.4, 1.5, 1.6, 2.0]
def o_filename(beta): return f"outputs/data_{beta}.csv"
def i_filename(beta): return f"inputs/params_{beta}.yml"
o_filenames = [o_filename(beta) for beta in betas]
processes = []
L = 100

measurements = 100
thermalization = 1000
record_rate = 100

refresh_rate = 0.1

# some functions

def parse_file(beta):
    data = []
    with open(o_filename(beta), 'r') as datafile:
        csv_data = csv.reader(datafile)
        for lines in csv_data:
            data.append([float(l) for l in lines])
    return np.array(data).T

class Observable(object):
    def __init__(self,label,func,scale = 'linear'):
        self.label = label
        self.func = func
        self.__uses_beta = len(signature(func).parameters) > 1
        self.scale = scale

    def __call__(self, data, beta):
        if self.__uses_beta:
            return self.func(data, beta)
        else:
            return self.func(data)

    def jackknife(self, data, beta, tau=0.5):
        cumsum = 0
        average = self(data, beta)
        for i in range(data.shape[1]):
            cumsum += ( self(np.delete(data,i,axis=1), beta) - average )**2
        return np.sqrt(2 * tau * cumsum)

global_observables = dict(
    # magnetization =     Observable(r"\langle|\bar\phi|\rangle", lambda data: np.mean(np.abs(data[0])) ),
    action          =   Observable(r"S", lambda data,beta: np.mean(data[0])/(beta*L**2) ),
    internal_energy =   Observable(r"E", lambda data,beta: np.mean(data[0])/(0.5*beta*L**2) ),
    chi_m           =   Observable(r"\chi_m", lambda data,beta: np.mean(data[1])/L**2, scale='log'),
    delta_m         =   Observable(r"\delta_m",
                            lambda data,beta:
                                2 * 10**5 * beta**4 * np.exp(-4*np.pi*beta) * np.mean(data[1])/(L**2),
                            scale = 'log')
)

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

@click.command()
@click.option('-q', '--quiet', help='run simulation quietly', is_flag=True)
def run(quiet):
    if not quiet:
        term = Terminal()
        if term.height+1<len(betas):
            raise Exception("terminal too short")
    for beta in betas:
        yaml_data = {
                'beta': beta,
                'L':    L,
                'measurements': measurements,
                'thermalization': thermalization,
                'record_rate': record_rate
                }
        with open(i_filename(beta), 'w') as yaml_file:
            yaml.dump(yaml_data, yaml_file)

        cmd = ["mpiexec",
                "-n", str(cores),
                "cpp_code/bin/sweep",
                "-i", i_filename(beta),
                "-o", o_filename(beta)]
        if quiet: cmd.append('-q')
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
    time.sleep(delay)


@click.command()
def view():
    chi_m_theory = [
            gvar(14.4, 0.4),
            gvar(23.5, 0.8),
            gvar(39, 3),
            gvar(82,8),
            gvar(191,21),
            gvar(370,48),
            gvar(871, 61),
            gvar(1634, 82),
            gvar(2285,58),
            gvar(2720,85)
        ]

    chi_m_betas = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]


    observ_names = ["internal_energy", "chi_m"]
    observ_objs = {o:global_observables[o] for o in observ_names}

    fig,axes = plt.subplots(len(observ_names),1,figsize=(12,6), sharex=True, squeeze=False)
    axes = [ax for row in axes for ax in row] # flatten axes

    axes[0].set_ylim((0,4.0))
    axes[0].set_title(f"L={L}, {measurements} measurements every {record_rate} sweeps, {thermalization} thermalization")

    axes[1].set_xlabel(r"$\beta$")

    plot_data = {k:[] for k in observ_objs.keys()}
    for beta, o_filename in zip(betas, o_filenames):
        data = parse_file(beta)
        for name,O in observ_objs.items():
            plot_data[name].append( gvar(O(data, beta), O.jackknife(data, beta)) )

    for ax, name in zip(axes, observ_names):
        d = plot_data[name]
        O = observ_objs[name]
        if name=="chi_m":
            ax.errorbar(chi_m_betas, [x.mean for x in chi_m_theory], yerr=[x.sdev for x in chi_m_theory], label="Berg & Lüscher", fmt='b.', capsize=3.)
        ax.errorbar(betas, [x.mean for x in d], yerr=[x.sdev for x in d], fmt='k.', label="Monte Carlo", capsize=3.)

        ax.set_ylabel('$'+O.label+'$')
        ax.set_yscale(O.scale)

        ax.set_xlim((0.4, 2.6))

    axes[1].legend()

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
def read_data(beta):
    fname = o_filename(beta)
    with open(o_filename, 'r') as datafile:
        csv_data = csv.reader(datafile)
        for lines in csv_data:
            data.append([float(l) for l in lines])
    return np.array(data).T

@click.command()
@click.option('-b', '--beta', required=True, type=float)
# @click.option('-O', '--observable', required=True, type=str)
def autocorrelation(beta):
    data = parse_file(beta)
    internal_energy =   data[0]/(0.5*beta*L**2)
    chi_m = data[1]/L**2

    name = f"Internal Energy: $L={L}$, $\\beta={beta}$"
    UWerr(internal_energy, name=name, plot=True, whole_plot=True)

@click.group()
def cli():
    pass

cli.add_command(run)
cli.add_command(view)
cli.add_command(autocorrelation)

if __name__ == "__main__":
    cli()
