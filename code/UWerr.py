'''
Python translation of UWerr.m from https://www.physik.hu-berlin.de/de/com/ALPHAsoft

Original author: Ulli Wolff
Transliterator: Stuart Thomas (snthomas01@email.wm.edu)
'''


import numpy as np
from math import sqrt, floor
from scipy.special import gammainc
from matplotlib import pyplot as plt



def UWerr(data, Stau=1.5, Nrep=None, name='NoName', quantity=0, plot=False, whole_plot=False):
    """

    UWerr(Data,Stau,Nrep,name,quantity,P1,P2,....)

    autocorrelation-analysis of MC time-series
    following the Gamma-method in
    ``Monte Carlo errors with less errors''
    by Ulli Wolff, hep-lat/0306017

    ----------------------------------------------------------------------------
     Ulli Wolff,   Nov. 2006 Version V6
     new: errobars on normalized autocorrelation (Luescher formula, see below)
          subleading terms in error of tauint changed again (see v4 of paper)
          protection against error^2 < 0 estimates
    ----------------------------------------------------------------------------
    input arguments beyond Data maybe omitted or set to []
    in the input list; they then get default value indicated in [D=..]


    Data     -- (N x Nalpha) matrix of measured (equilibrium!) data
                 N = total number of measurements
                 Nalpha = number of (primary) observables
                 if your data are in a different format, consider
                 the numpy commands `reshape' and `permute'

    Stau     -- guess for the ratio S of tau/tauint [D=1.5]
                if set = 0, absence of autocorrelations is *assumed*

    Nrep     -- vector [N1 N2 N3 ...] specifying a breakup of the N rows
                of Data into replica of length N1,N2,.. (N=sum(Nrep)!)  [D=[N]]
                The replica distribution is histogrammed and a Q-value
                (=probability to find this much or more scatter) is given
                for R >= 2 replica; if you have one history you may set Nrep to
                artificially split it into >=2 replica to see their distribution

    name     -- if string: name of observable for titles of generated plots
                if not string: all plots are supressed [D='NoName']

    quantity -- either:
                scalar function for the derived
                observable F; it has to operate on a row-vector of length
                Nalpha as first argument; optional parameters P1,P2,... are
                passed on to this function as 2nd, 3rd, .... argument
             -- or:
                integer between 0 and Nalpha-1 to analyze primary observable [D=0]

    ----------------------------------------------------------------------------
    output::
     value   -- estimate for F
     dvalue  -- its error (incl. autocorrelation effects)
     ddvalue -- statistical error of the error
     tauint  -- integrated autocorrelation time
     dtauint -- statistical error of tauint
     Qval    -- Q-value of the replica distribution if R >= 2
                (goodness of fit to a constant)
    ----------------------------------------------------------------------------
    """

    if data.ndim==1:
        data  = np.expand_dims(data, axis=1)
    N,Nalpha = data.shape

    if Nrep is None:
        Nrep = np.array([N])


    if np.any(Nrep != np.round(Nrep)) or np.any(Nrep < 1) or np.sum(Nrep) != N:
        raise Exception("inconsistent N,Nrep")
    if not callable(quantity):
        if quantity != round(quantity) or quantity < 0 or quantity >= Nalpha:
            raise Exception("illegal numerical value of quantity")
        primary = True
        primaryindex = quantity
    else:
        primary = False

    R = len(Nrep) # number of replica

    # means of primaries
    abb = np.mean(data, axis=0) # total mean, primary obs.
    abr = np.zeros((R, Nalpha)) # replicum-mean, primary obs.
    i0=1
    for r in range(R):
        i1 = i0 - 1 + Nrep[r]
        abr[r,:] = np.mean(data[i0:i1,:], axis=0)
        i0=i0+Nrep[r]

    if primary:
        Fbb = abb[primaryindex]
        Fbr = abr[:, primaryindex]
    else:
        Fbb = quantity(abb)  # total mean, derived obs.
        Fbr = zeros((R,1))   # replica-means, derived obs.
        for i in range(R):
            Fbr[i]=quantity(abr[i,:])

    Fb = np.sum(Fbr * Nrep) / N # weighted mean of replica means

    # for the gradient of f and project fluctuations
    if primary:
        delpro = data[:,primaryindex] - abb[primaryindex]
    else:
        fgrad = np.zeros(Nalpha)
        h = np.std(data)/sqrt(N)
        ainc=abb
        for alpha in range(Nalpha):
            if h[alpha] == 0:
                fgrad[alpha]=0
            else:
                ainc[alpha] = abb[alpha] + h[alpha]
                fgrad[alpha] = quantity(ainc)
                ainc[alpha] = abb[alpha] - h[alpha]
                fgrad[alpha] = fgrad[alpha]-quantity(ainc)
                ainc[alpha] = abb[alpha]
                fgrad[alpha] = fgrad[alpha] = fgrad[alpha]/(2*h[alpha])


        # projected deviations
        delpro = data @ fgrad.T - np.dot(abb, fgrad)

    value = Fbb
    dvalue = 0
    ddvalue = 0
    tauint = 0.5
    dtauint = 0
    Qval = np.array([])


    # compute Gamma, automatic windowing

    if Stau == 0:
        Wopt = 0
        tmax = 0
        flag = False
    else:
        tmax = floor(min(Nrep)/2)
        flag = True
        Gint = 0
    # values of W=0

    GammaFbb = np.empty(tmax+1)
    GammaFbb[0] = np.mean(delpro**2)
    # sick case:
    if GammaFbb[0] == 0:
        print("WARNING: no fluctations")
        return value, dvalue, ddvalue, tauint, dtauint, Qval


    for t in range(1, tmax+1):
        GammaFbb[t] = 0
        i0 = 0
        for r in range(R):
            i1 = i0 - 1 + Nrep[r]
            GammaFbb[t] = GammaFbb[t] + np.sum(delpro[i0:i1-t] * delpro[i0+t:i1])
            i0=i0+Nrep[r]


        GammaFbb[t] = GammaFbb[t]/(N-R*t)
        if flag:
            Gint=Gint+GammaFbb[t]/GammaFbb[0]
            if Gint <= 0:
                tauW= np.finfo(float).eps
            else:
                tauW = Stau/np.log((Gint+1)/Gint)
            gW = np.exp(-t/tauW) - tauW/sqrt(t*N)
            if gW < 0: # this w is taken as optimal
                Wopt = t
                tmax = min(tmax, 2*t)
                flag = False
    if flag:
        print(f"WARNING: windowing condition failed up to W = {tmax}")
        Wopt=tmax

    CFbbopt = GammaFbb[0] + 2*np.sum(GammaFbb[1:Wopt])
    if CFbbopt <=0:
        print(f"WARNING: Gamma pathological: estimated error^2 < 0")
        return value, dvalue, ddvalue, tauint, dtauint, Qval

    GammaFbb = GammaFbb+CFbbopt/N
    CFbbopt=GammaFbb[0] + 2*np.sum(GammaFbb[1:Wopt])
    sigmaF=np.sqrt(CFbbopt/N)
    rho=GammaFbb/GammaFbb[0]
    tauintFbb=np.cumsum(rho)-0.5

    # bias cancellation for the mean value

    if R >= 2:
        bF = (Fb-Fbb)/(R-1)
        Fbb=Fbb - bF
        if abs(bF) > sigmaF/4:
            print(f"WARNING: Gamma pathological: estimated error^2 < 0")

        Fbr = Fbr - bF * N / Nrep
        Fb = Fb - bF * R




    # answers to be returned

    value = Fbb
    dvalue = sigmaF
    ddvalue = dvalue*sqrt((Wopt+0.5)/N)
    tauint = tauintFbb[Wopt]
    dtauint = tauint * 2 * sqrt((Wopt-tauint+0.5)/N)


    # Q-value for replica distribution if R >= 0

    if R >= 2:
        chisq = np.sum((Fbr - Fb)**2 * Nrep)/CFbbopt
        Qval = 1 - gammainc((R-1)/2, chisq/2)
    else:
        Qval = np.array([])

    # plotting
    if plot:
        if Stau !=0:
            # plot of GammaF(t)/Gamma(0)=rho(t)
            # contruct errors acc. to hep-lat/0409106 eq. (E.11)
            # pad zeros to simplify summation:x
            rho[tmax+1:]=0

            drho = np.empty(tmax+1)
            for t in range(tmax+1):
                k = np.arange(max(0, t-Wopt-1), t+Wopt-1)
                drho[t] = sum((rho[k+t] + rho[abs(k-t)] - 2*rho[t]*rho[k])**2)
                drho[t] = sqrt(drho[t]/N)

            fig, (rho_ax, tau_ax) = plt.subplots(2,1, figsize=(16,14))
            rho_ax.axhline(0., ls='--', c='r')
            rho_ax.axvline(Wopt, c='r')
            rho_ax.plot(np.arange(tmax+1), rho[:tmax+1])
            rho_ax.fill_between(np.arange(tmax+1), rho[:tmax+1]-drho, rho[:tmax+1]+drho, alpha=0.3, color='g')
            rho_ax.set_title('normalized autocorrelation of ' + name)
            rho_ax.set_ylabel(r"$\rho$")
            rho_ax.set_xlim((0, tmax))

            tau_ax.axhline(tauint, ls='--', c='r')
            tau_ax.axvline(Wopt, c='r')
            tau_ax.plot(np.arange(tmax+1), tauintFbb[:tmax+1])
            dtau = tauintFbb[:tmax+1] * np.sqrt(np.arange(tmax+1)/N)*2
            tau_ax.fill_between(np.arange(tmax+1), tauintFbb[:tmax+1]-dtau, tauintFbb[:tmax+1]+dtau, alpha=0.3, color='g')
            tau_ax.set_title('$\\tau_{int}$ with statistical erros of ' + name)
            tau_ax.set_ylabel(r"$\tau_{int}$")
            tau_ax.set_xlabel(r"$W$")
            tau_ax.set_ylim((0.5, None))
            tau_ax.set_xlim((0, tmax))
            plt.show()






    return value, dvalue, ddvalue, tauint, dtauint, Qval

