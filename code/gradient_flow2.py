from parallel_utils import RANK
from mcmc.keys import WOLFF, SWENDSEN_WANG
from mcmc import Phi4Lattice, RandomWalk, GradientFlow, Recorder
from mcmc.recorder import binder_cumulant
import numpy as np
from time import time
import sys
if 'view' in sys.argv:
    from matplotlib import pyplot as plt
import pickle
import gvar
from datetime import timedelta
from functools import partial
import os
SEED = True
Ls = [int(sys.argv[1])]
lam = 0.5
m02s = [-0.80,-0.72, -0.64]
taus = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]

# sweeps = 10**5
cluster_method = WOLFF
thermalization = 10**4
record_rate = 100
cluster_rate = 5

measurements = 10**3

sweeps = thermalization + record_rate * measurements
if RANK==0:
    print(f"{sweeps} sweeps")




DATAPATH = f'data/{os.path.basename(__file__)}.pickle'

def main():

    GF = GradientFlow()
    all_recorders = []
    for m02 in m02s:
        if RANK==0: print(f"m02 = {m02:.2f}")
        for L in Ls:

            if RANK==0:
                hooks = [partial(GF.flow_evolution, tau=tau, action=True) for tau in taus]
                recorders = [Recorder(thermalization=thermalization, rate=record_rate, gif=False) for tau in taus]
            else:
                hooks = []
                recorders = []


            l = Phi4Lattice(dim=L, m02=m02, lam=lam)
            rw = RandomWalk(l)


            if RANK==0:
                start = time()

            rw.run(sweeps, cluster_method=cluster_method,  cluster_rate=cluster_rate, hooks=hooks, recorders=recorders, progress=True)

            if RANK==0:
                dt = time() - start
                print(f"Executed in {int(dt//60)}:{str(int(dt%60)).zfill(2)}")
                # if recorder.gp is not None: recorder.save_gif("plots/BC.gif")
            all_recorders.extend(recorders)

    if RANK==0:
        pickle.dump(all_recorders, open(DATAPATH, 'wb'))
        print(f"Saved pickle to {DATAPATH}")





def view():
    recorders = pickle.load(open(DATAPATH, 'rb'))

    quantities = ['magnetization','susceptibility','binder_cumulant', 'action']

    # for i,r in enumerate(recorders):
        # n = r.record_count
        # assert n==101
        # phi2_avg = sum(r.values['phi2'])/n
        # phi4_avg = sum(r.values['phi4'])/n
        # plt.plot(r.values['phi4'], label='$\phi^4$')
        # plt.plot(r.values['phi2'], label='$\phi^2$'):w
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

    fig, axes = plt.subplots(5,1, figsize=(16,10))
    L = Ls[0]

    for i,m02 in enumerate(m02s):
        for ax, quantity in zip(axes[:-1], quantities):
            some_recorders = recorders[i*len(taus):(i+1)*len(taus)]

            # test_r = some_recorders[-1]
            # phis = np.array(test_r.values['phi'])

            # gvs = test_r.gvars
            # phi4 = gvs['phi4']
            # phi2 = gvs['phi2']
            # plt.figure()
            # plt.plot(phis**2, label='$\phi^2$')
            # plt.plot(phis**4, label='$\phi^4$')
            # plt.legend()
            # plt.show()

            # print( 1 + phi2.sdev**2 / phi2.mean**2)
            # print( phi4.mean / phi2.mean**2)
            # print(gvar.evalcov([phi4, phi2**2]))
            # print(np.sqrt(gvar.evalcov([phi4, phi2**2])))

            # phi2_ = gvar.gvar(phi2.mean, phi2.sdev)
            # phi4_ = gvar.gvar(phi4.mean, phi4.sdev)

            # print(phi4, phi2**2)
            # print(phi4/ phi2**2)
            # print(phi4/ phi2**2)
            # print(phi4_/ phi2_**2)
            # print(1 - gvs['phi4'] / ( 3 * gvs['phi2'] ** 2))
            # quit()

            means = np.array([r.derived_values[quantity] for r in some_recorders])
            stds  = np.array([r.derived_errors[quantity] for r in some_recorders])

            ax.errorbar(taus, means, yerr=stds, fmt='o', label=f"$m_0^2={m02}$", capsize=5,zorder=i)
            ax.legend()
            if quantity=="binder_cumulant":
                ax.axhline(2/3, c='r')
                ax.axhline(0, c='k')
            # if quantity=="magnetization":
                # ax.set_ylim( (0.0, np.sqrt(-m02s[-1]/lam)) )
            # print("phi2",[r.values['phi2'] for r in some_recorders]) # EDIT

            ax.set_ylabel(Recorder.derived_observables[quantity].label)



       #  bimod_axis = axes[-1]
        # bimod_axis.plot(m02s, [bimodality(r.values['phi']) for r in some_recorders], 'o', label=f"$N={L}$")
        # bimod_axis.set_ylabel("B")


        # time_ax = axes[-1]
        # some_tds = tds[i::len(Ls)]
        # time_ax.plot(m02s, some_tds, 'o', label=f"$N={L}$")
        # time_ax.set_ylabel("Execution time (m)")
        # time_ax.set_yscale("log")
        # time_ax.grid(b=True, which='both')
        # time_ax.set_ylim(10**3, 10**6)
    tds = [0.]
    execution_time_str = f"{sum(tds)//60:.0f}h{sum(tds)%60:.0f}m"
    axes[0].set_title(f"{measurements} measurements every {record_rate} sweeps, {thermalization} thermalization, Jackknife errors, (execution time: {execution_time_str})")

    plt.show()

if __name__=="__main__":

    if 'view' in sys.argv:
        view()
    else:
        main()
