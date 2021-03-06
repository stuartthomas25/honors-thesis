from parallel_utils import RANK
from mcmc.keys import WOLFF, SWENDSEN_WANG
from mcmc import Phi4Lattice, RandomWalk, Recorder
from mcmc.recorder import binder_cumulant
import numpy as np
from time import time
import sys
if 'view' in sys.argv:
    from matplotlib import pyplot as plt
import pickle
from datetime import timedelta

SEED = True
Ls = [64]
lam = 0.5
m02s = [-0.80,-0.76, -0.72, -0.68, -0.64]
# sweeps = 10**5
cluster_method = WOLFF
thermalization = 10**3
record_rate = 100
cluster_rate = 5

measurements = 10**3

sweeps = thermalization + record_rate * measurements
if RANK==0:
    print(f"{sweeps} sweeps")

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



def main():

    all_recorders = []
    for m02 in m02s:
        if RANK==0: print(f"m02 = {m02:.2f}")
        for L in Ls:

            if RANK==0:
                recorders = [Recorder(thermalization=thermalization, rate=record_rate, gif=False)]
            else:
                recorders = []

            l = Phi4Lattice(dim=L, m02=m02, lam=lam)
            rw = RandomWalk(l)


            if RANK==0:
                start = time()

            rw.run(sweeps, cluster_method=cluster_method,  cluster_rate=cluster_rate, recorders=recorders, progress=True)

            if RANK==0:
                dt = time() - start
                print(f"Executed in {int(dt//60)}:{str(int(dt%60)).zfill(2)}")
                # if recorder.gp is not None: recorder.save_gif("plots/BC.gif")

            all_recorders.extend(recorders)

    if RANK==0:
        pickle.dump(all_recorders, open('data/binder_cumulant.pickle', 'wb'))

        print("Saved pickle")




def view():
    recorders = pickle.load(open('data/binder_cumulant.pickle', 'rb'))

    quantities = ['magnetization','susceptibility','binder_cumulant']

    # for i,r in enumerate(recorders):
        # n = r.record_count
        # assert n==101
        # phi2_avg = sum(r.values['phi2'])/n
        # phi4_avg = sum(r.values['phi4'])/n
        # plt.plot(r.values['phi4'], label='$\phi^4$')
        # plt.plot(r.values['phi2'], label='$\phi^2$')
        # plt.legend()
        # plt.show()
        # print(phi2_avg**2, phi4_avg)
        # print(i, 1 - phi4_avg/(3*phi2_avg**2))
        # print()

    # fig, axes = plt.subplots(4,1, figsize=(16,10))
    # for ax, quantity in zip(axes, Recorder.primary_observables):
        # for i,L in enumerate(Ls):
            # recorder = recorders[i+15]
            # values = recorder.values[quantity]
            # ax.plot(values, label=f"$N={L}$")
            # # ax.errorbar(m02s, derived_values, yerr=stds, label=f"$N={L}$")
            # ax.legend()

            # ax.set_ylabel(quantity)
    # plt.show()
   #  for r in recorders:
        # phi = r.gvars['phi']
        # phi2 = r.gvars['phi2']
        # phi4 = r.gvars['phi4']
        # print(1,np.sqrt(phi2.mean - phi.mean**2), phi.sdev)
        # print(2,np.sqrt(phi4.mean - phi2.mean**2), phi2.sdev)

    print(len(recorders))
    for r in recorders:
        r.finalize_values()

    tds = []
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

    fig, axes = plt.subplots(4,1, figsize=(16,10))
    for i,L in enumerate(Ls):
        for ax, f in zip(axes[:-1], [magnetization, susceptibility, binder_cumulant]):
            some_recorders = recorders[i::len(Ls)]
            quantity = f.__name__

            means = []
            stds = []
            for r in some_recorders:
                phi = r.values['phi']
                means.append(f(phi))
                stds.append(jackknife(f, phi))

            # means = np.array([r.derived_values[quantity] for r in some_recorders])
            # stds  = np.array([r.errors[quantity] for r in some_recorders])
            # stds  = np.array([r.derived_errors[quantity] for r in some_recorders])

            ax.errorbar(m02s, means, yerr=stds, fmt='o', label=f"$N={L}$", capsize=5,zorder=i)
            # ax.legend()
            if quantity=="binder_cumulant":
                ax.axhline(2/3, c='r')
                ax.axhline(0, c='k')
            if quantity=="magnetization":
                ax.set_ylim( (0.0, np.sqrt(-m02s[-1]/lam)) )
            # print("phi2",[r.values['phi2'] for r in some_recorders]) # EDIT

            ax.set_ylabel(Recorder.derived_observables[quantity].label)



        bimod_axis = axes[-1]
        bimod_axis.plot(m02s, [bimodality(r.values['phi']) for r in some_recorders], 'o', label=f"$N={L}$")
        bimod_axis.set_ylabel("B")

    execution_time_str = f"{sum(tds)//60:.0f}h{sum(tds)%60:.0f}m"
    # axes[0].set_title(f"{measurements} measurements every {record_rate} sweeps, {thermalization} thermalization, Jackknife errors, (execution time: {execution_time_str})")
    print(f"{measurements} measurements every {record_rate} sweeps, {thermalization} thermalization, Jackknife errors, (execution time: {execution_time_str})")
    axes[-1].set_xlabel(r"$m_0^2$")

    plt.show()

if __name__=="__main__":

    if 'view' in sys.argv:
        view()
    else:
        main()
