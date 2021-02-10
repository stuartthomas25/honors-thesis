---
title: 'Thesis Notes'
author:
    - Stuart Thomas
    - Chris Monahan
bibliography:
    - zotero.yaml
link-citations: true
nocite : |
    @*
linestretch: 1.5
...

## Monday, June 15
-  Preliminary meeting
-  We will begin with $\phi^4$ term due to lower energy bound.
-  Beginning with code in Python, switch to C/C++ if necessary.

## Tuesday, June 16
-  Preliminary concepts to understand:
    -  scalar field theory on lattice
    -  Markov chains and Monte Carlo
    -  the gradient flow
    -  $O(n)$ symmetry
-  General components to research, executed in parallel
    -  Reading
    -  Mathematical analysis
    -  Writing code
    -  Present (writing, plot generation,…)
-  Things to get out of dissertation [@schaich2006]
    -  Markov chain
    -  Cluster algorithms
        -  Metropolis is a method for narrowing possibilities by
            accepting only some changes. It can get stuck in a local
            mininum, loss of ergodicity. We solve these with cluster
            algoriths. Wolff grows clusters probabilistically and flips,
            while Swenson and Yang identifies clusters and flips them
            probabilistically.
        -  Near a phase transition, correlation length grows and
            changes become less likely to be accepted: need clusters.
            Clusters dont work far from the phase transition. This is
            manifested as a sequence of a few metropolis steps and a
            cluster step.
-  Research plan
    -  Start with 2D $\phi^4$
        -  Set up lattice with sign flip for reflection
        -  Use Markov chain Monte Carlo to simulate.
        -  Measure, magnetization and suseptibility, Binder cumulant
    -  Transition to 3D, then maybe transition to C/C++.
    -  Implement the gradient flow
    -  move to 2-3d nonlinear sigma model.
    -  Motivation: the nonlinear sigma model works for QCD given the
        asymtotic freedom. We may also want to explore topology.

## Wednesday, June 17
-  Code tips:
    -  Try Swendsen-Wang algorithm in addition to Wolff
    -  Print out time taken
    -  Optimize Hamiltonian
    -  Save every tenth measurement or store configurations to
        calculate path integral. Exclude thermaliztion (first 200)
    -  Write in terms of sweeps, not iterations
    -  Parallelize (look into checkerboard algorithm)
    -  Implement Binder cumulant and suseptibility
    -  Store every few states
    -  Look into multigrid algorithm
-  Reading on Monte Carlo Markov chain and cluster algorithms.

# Friday, June 19
-  Looking over code
    -  Might be too slow, move to C/C++ eventually?
    -  Why is the energy increasing with metropolis algorithm?
    -  Shift to using action instead of Hamiltonian.
    -  Transition from Broken Phase, look at $\mu_0^2$ term.
    -  Profile code for possible optimizations.

## Monday, June 22
-  Coding
    -  Add plots to `.gitignore`.
    -  Try parallelizing code, using either multigrid or checkerboard
        algorithm
-  Reading
    -  Start to focus more on understanding the theory behind research.
    -  Read dissertation [@schaich2006] Chap. 6.5 and 6.6, take notes

        on questions.
    -  Newman [@newman1999] (Main textbook for Monte Carlo in

        Statistical Physics)
-  Just some things to remember
    -  Correlation functions correlate values in statistical systems
        and relate to propagators in QFT.
    -  The problem of renormalization: $a\rightarrow0$ leads to
        unbounded correlation function. As real physical lengths are
        measured in terms of the lattice constant, these sometimes tend
        to infinity.

## Wednesday, June 24
-  Code
    -  Parallelize
    -  Transition to numpy
-  Theory Question
    -  renormalize and regularlize: what do they mean?
    -  Look at LePage

## Tuesday, June 30
-  Code
    -  Try `mpi4py`.
    -  Move to 3D (this may decrease parallel overhead)
-  Reading
    -  Continue Reading Collins, others.

## Thursday, July 2
-  Coding
    -  Continue implementing MPI
    -  Parallelization may be more apparent in 3D
-  Reading
    -  Dirac fermions, represented by 4D spinor field (spin up/down,
        electron/positron)
    -  Look at Tong (Chap. 4)
    -  Charge (See Tong, Noether’s Theorem)
-  Next week back: start gradient flow on linear $\phi^4$ model

## Monday, July 13
- Update from last week
    - Implemented MPI, slowdown may be due to thermal throttling?
    - Questions:
        - Noether's Theorem, 4-current?
        - Star vs. dagger

- Invariance vs covariance
    - $Q$ (charge) is invariant, not covariant. Derivative is 0 (conserved) so no effect of boost. See (StackExchange)[https://physics.stackexchange.com/questions/270296/what-is-the-difference-between-lorentz-invariant-and-lorentz-covariant]

- Operator product expansion    
    - Taylor (Laurent in reality) series for operators
    - Used to expand nonlocal (slightly) oeprators using local operators



- TODO:
    - Read conference proceedings, then paper with Orginos
    - Define $\rho$ field (Eq. 2.4), implement it using the exact solution (Eq. 2.5). This will require a FFT

## Wednesday, July 15
- Code organization
    - Perhaps implement the nonlinear $\sigma$ model (at a later date) and the $\phi^4$ model as subclasses of a generic lattice class with abstract methods of action, etc.

- Gradient flow
    - Splitting momentum at half of lattice is a-ok.
    - Next, implement observable!

- New terms:
    - autocorrelation times (see Schaich -> MC Textbook -> Wolff): not a QFT concept!
    - Gamma analysis: used to calculate autocorrelation
    - summation window: part of the Gamma analysis

- Three classes of QFTs
    - renormalizable: infinities can be absorbed by a finite number of counterterms
    - nonrenormalizeable: requires infinitely many counterterms
    - super-renormalizeable: there is only one parameter that is divergent

## Friday, July 17

- Critical mass problem
    - Check multiple measures (Binder Cumulant, bimodality, etc.)
    - https://journals.aps.org/prd/abstract/10.1103/PhysRevD.58.076003
    - Infinite volume limit?
    - Turn off coupling?
    - Behavior for one node? Two nodes?
    - Metropolis checkerboard? (see [@schaich2006, pg. 79])
    - Look at Schaich's calculations
    - Note: No true phase transition in finite system since phase transitions are defined by correlation lengths going to infinity

## Monday, July 20

- Remarks on critical point graph
    - I should be averaging the measurements from run.
    - Probably will need many more measurements.
    - Spend more time thermalizing
    - Calculate autocorrelation time, use to determine the `record_rate`
    - Production quality: try 10,000 sweeps, with 1000 thermalization
    - Run with slightly larger lattice to verify critical point

- Todo:
    - Flow time dependence (with action)
    - Average all measurements (incorporate this into `recorder.py`)
        - Plot mean and standard error as error bars
    - Autocorrelation times using Gamma analysis
        - Paper has Mathematica notebook link, though it may be broken.
    - Double check that large lattice volumes show correct critical mass.

## Wednesday, July 22

- Things I did:
    - Reorganized the `recorder.py` to include measurements and means/errors.
    -

- Questions:
    - Implementing the gradient flow: how does this play into the Monte Carlo simulation?


- TODO:
    - Critical mass
        - Plot histogram for two points of single lattice (thermalized)
        - See if this is affected by a large thermalization, cold start
        - Binder cumulant should be step function

    - Autocorrelation times
        - See Eq. 41
        - Get MatLab code
    
    - Gradient flow
        - Each measured lattice should be evolved in flow time before recording

## Friday, July 24
 
- Things to do:
    - Gradient flow
    - Try $m_0^2 = -0.8$, see if it become bimodal with large thermalization. 
    - Move towards critical value: does $\tau_{int}$ increase? (This is what we expect)
    - Perhaps change recording rate?

- Is it time to transcibe the code?
    - Explore efficiency boost
    - C++?
    - Cython?
    - Local cluster?



## Monday, July 27

- Shift to utilizing the 0th node for computation.

- Things to note:
    - Histogram in conference proceedings bins single lattice, not average over Markov chain.

- Questions:
    - collecting plots

- To-Do:
    - Run a proper histogram for different thermalizaitons, see if 10^4 is necessary.

## Wednesday, July 29

- Plots look good, but are a little wide and not localized on the convergent point
    - Try a cold start
    - Look at 256
    - Average a couple trials

- To-Do:
    - See if I can find the width issue
    - Check autocorrelation times (with/without cluster)
        - close to critical value
    - wait on gradient flow until we are more confident in Markov chain.


## Friday, July 31

- Calculate autocorrelation times (as a function of mass), binder cumulant
- Check Wolff algorithm, try Swendsen Wang
- Multigrid algorithm?

- Focus on metropolis algorithm. We don't want to study a cluster algorithm on the backdrop of a broken metropolis algorithm.

- Literature search: must the metropolis algorithm have randomly included sites?
    

## Thursday, August 6

- Progress:
    - Autocorrelation: see page 6 from Schaich
        - Affected by the measurement rate?
    - Wolff implemented in C, large speed-up
    - Fixed neighbor bug in Wolff algorithm
    - Histogram issue still persists


## Friday, August 7

- Binder cumulant shows very different critical mass

- Going forward:
    - Start taking personal notes on the research with sketches, ideas, progress
    - Always calculate all of the observables


- Try the Binder Cumulant without Wolff
- See if the result depends on thermalization
    - This result should not be the same as 0 thermalization.



## Wednesdesday, August 12

- There are some major problems with my averaging.
- TODO: implement this fix

Graduate School:

- Aim high
- Best possible school that you're happy at
- Name recognition matters
- If you can afford to do the applications, apply to a lot
- Maybe not scientific instrument lab?

- You will have a choice betweeen famous advisors:
    - Famous: letters and papers will carry weight
    - Not-famous: much less attanetion


- Are you emailing your future advisor?
    - depends on place
    - At MIT, less important to single out a specific professor.
    - At a small school (W&M), you need to identify someone before hand.
        - People can actually change advisors
- Possible schools:
    - UCSD
    - Berkley
    - Rutgers
    - UMD
    - UCSB
    - MIT
    - Yale


## Friday, August 14
- First, check previous version and make sure the numbers are not the same.
- Try saving some lattices, then calculating observables in Mathematica

## Tuesday, August 25
- Progress:
    - Fixed BC issue by fixing the order of calculations. New order:    
        1. calculate volume-average
        2. take any exponents
        3. take the ensemble average
        4. compute any derived quantities

    - This is supported by Eq. 4.19 [@newman1999].

- Future work

    - Calculate errors for derived quantities using `gvar`.
    - Transition to GitHub (ID is `cjmonahan`).
    - Run through measurements that you can think of. Try other measurements (e.g. bimodality) and generate definitive list of plots.
    - Create plots without connected lines


## Thursday, August 27

- Problem with Binder-Cumulant: errors are too large.
    - Implement Jackknife method (see [@toussaint1989]) to measure statistics of variables.
        - Try this on 128 lattice, maybe once overnight
    - Migrate to cluster
    - Or perhaps beforehand, try gradient flow.


## Tuesday, September 1
- Questions:
    - Do I apply it to every measurement?
        - Yes
    - Are flow time and Monte Carlo time independent?
        - Not really. This is actually where the gradient flow originates, but its not the point of our research.
    - What is the gradient flow scale?
        - We know that $mass \times length$ is dimensionless, so we can measure the flow time in $mass^{-1}$.
        - In other more realistic theories, we can write the flow time as something physical, like the proton size.
        - In scalar field theory, we can do the same thing, but there is no physical analog.
        - The length scale of flow time should be written in terms of $\lambda$ or $m_0^2$, with $m_0^2$ being the slightly more natural choice.

- Todo:
    - Pick three fixed masses, study the flow time at each.
    - Evolve in flow time for each measurement.
    - Try this at $N=96$.
    - Add the action as an observable.


## Friday, September 4
- Gradient flow should have a 0 imaginary component, so I can take the real component safely, though I should check this.


- Todo:
    - Add momentum constant to Gradient Flow evolution function.
    - Fix flow evolution bug!

## Tuesday, September 8
- Issue: Observable values are actually completely flat, may be a memory error in the code. This needs debugging.
- Plot the results on a smaller flow time scale.


## Friday, September 11
- Reminder: Must rescale momenta by $\pi/L$.
- Problem: Observables look flat

- We won't reproduce Figure 3 from proceedings.
- Continue to nonlinear sigma model.
    - Use [@aoki2015] to find other sources
    - Use inspirehep to find other sources.
    - Check out [http://www.scholarpedia.org/article/Nonlinear_Sigma_model]

## Friday, September 18

- For Wolff algorithm, chose arbitrary vector to flip along.

## Tuesday, September 22

- Converting Cython to C++
    - Converting to flat lattice
    - Swendsen & Wang? we should implement this if Wolff works.

- The nonlinear $\sigma$ lagrangian:
    - $g$ value is related to $\sqrt{\lambda}$.
    - also does affect system since factor of action is significant in path integral formalism.
    - No mass term.

## Friday, September 25


- 
- Questions:
    - New phi value, guarantee magnitude less than 1 with the other 3 components?
        - Should be using rotation matrices.
        - Take a look at SU(2) Metropolis literature and see how they generate norm-preserving rotation matrices computationally.
    - What exactly are we calculating?
        - Twist-2 operators
        - Condensates
            - VEV in non-perturbative theory.
        - Topology?
        - These could all be topics of our research
    - Prof. Chris: write notes next week about why the nonlinear $\sigma$ model is "cool".


\newpage

## Wednesday, September 30

- Implement MPI in C++.

- There actually is an application of nonlinear-sigma, in fact multiple. 
    - Heisenbuerg Ferro Magnetic, real application to condensed matter.
    - Applications to string theory.
    - Same properties as Yang-Mills Gauge Theories
        - O(2) renormalizable, nontopological solutions of classical EOM at finite action
        - Asymtoticaly free
        - Mass gap in nonperturbative theory
        - Large N limit
        - Dynamical generation of bosons, (like phonons, not explicit in Lagrangian).
    - In $d=4$, the NLSM is starting point for chiral perturbation theory.

## Tuesday, October 6
- Added MPI to C++
    - Some issues with conversion to C array

- Make sure in implementation to keep lattice data continuous.

## Wednesday, October 14
- How to debug compilation errors?
    - Can't use print statements like Python.
    - I can use dummy routines to isolate issues.

- Priorities:
    - Ensure that Enrico's project is on the Arxiv while still doing work on my thesis in good faith.

## Wednesday, October 21

- Prof. Chris will ask Prof. Orginos about static memory allocation for rank and size.
- Code finally compiles!
- Action looks good but magnetization is incorrect.
    - Get other measurements (suscpetibility, BC, etc.)
    - No need to export lattice average of $\phi^2$.

## Wednesday, October 28
- Progress
    - Plot with $U$ and $\chi$
    - I think I can statically allocate data in the Sweeper class using an initializer list.

- Do we need MPI?
    - Yes, in 3D it may be necessary since we won't have as many points.
    - Also is necessary for psychological help
    - Also allows adaptive sampling

- Mid Year report due Monday



## November 4 

- Comparing Python and C++ part. Some issues:
    - Susceptibilities have errors that are too large.
    - Ratio of susceptibilities also wrong, expect a sharper dropoff.
    - Binder Cumulant has wrong final values, susceptibilities.
    - A good future solution is plotting the difference and a horizontal 0 line, also plot error.
    - I should rerun everything with the same data.
    - Also put the data in a nice plot.
    - Make note of the speed up.


- Transition to $O(3)$ model.
- I realized that vectors are continuous memory, almost no slowdown.

## November 11

- Updates:
    - Ran calculations again
        - Execution time is actually on a similar order or magnitude
        - MPI on 4 cores actually slows progress
        - A couple of outliers, checked closer, seems right
    - Redeveloped gif debugging tool for C++
        - Shows similar lattice for both methods.
        - Can be ported to $O(3)$ model
- Todo:
    - For difference plot, use larger number of measurements.
    - Still discrepancy with susceptibility error bars
        - See if this changes with lattice size, ensemble size.
        - Jackknife methods
        - Export $\bar\phi$ values in python code, compare variances with independent code.

    - Presentation slides by Friday
    - Email Averett about a final draft for the mid-year report.
    - Make changes to the mid-year report 

---

# Spring Semester

## January 12

- Some clarifying theory questions:
    - Twist-2 operators: not trivial! mass dimension is not just mass of the field. 

- __Shift project towards topology__
    - More interesting
    - More applications to condensed matter
    - See references [@bietenholz2018] and [@mejia-diaz2018].
        - Unfortunately, these do exactly what we were planning to do, however they do demonstrate that this research is paper-worthy.

- Todo in the meantime:
    - clean up code

- Writing report:
    - Try to start early
    - Specifically, fix "long-term" issues in mid-year report.


## January 21


- New direction, since papers cover what we did:
    - Look at instantons?
- __Topology__
    - What is Q?
        - Global operator
        - $\chi$: topological suscpetibility

    - Wedge operator ($\wedge$).
        - Does it act in real space or spin space?
        - Must do further research.


- **General research question**: How do I go about figuring out mathematical expressions.
    - Keep Googling phrases! 
        - try the arxiv!
    - Find other papers with the same expression.
    - Compare different definitions from different papers.
    - Ask someone.

- How to determine gradient flow in $O(3)$ model?
    - Runge-Kutta: check out Numerical Recipe
    - Try $O(4)$ approx. 

- Todo:
    - First order: gradient flow and topological charge
    - Second order: susceptibilty


# January 27
- Parallel Runge-Kutta method
    - Prof. Chris will ask Prof. Orginos
    - Go ahead and implement serial algorithm, then check if adaptive step size changes efficiency significantly.
    - Implement second-order Runge-Kutta and compare efficiency/accuracy.
    - Matrix options: look at parallelism across the method.

- What statistics to use
    - Maybe use the term "observable"
    - Magnetization may be important, as with action and susceptibility
    - Topological density
    - For topological charge see [@berg1981]


- Todo:
    1. Runge-Kutta gradient flow
    2. Adaptive step size
    3. Check internal energy from [@berg1981]
    4. Then, consider topological charge/suscpetibility (Fig 2, [@berg1981])

#  February 1
- Reproduction of [@berg1981]:

![Internal Energy](internal_energy)

- Gradient flow:
    - Action increases with flow time
    - May be due to normalization.

- Orginos: Runge-Kutta parallel python code.
    May be a shortcut to calculate the gradient flow :[@luscher2010]:
        - Used with $SU(3)$ but not $O(3)$.
        - Could be mapped but the mathematics are very difficult (see general technique in [@munthe-kaas1999]).
        - Note that Bietenholz does not use L\"uscher's method.

- Internal Energy:
    - Run 100x100 and compare with specific $\beta$ values from Table 1 [@berg1981].
    - Try some more $\beta$ values, including 2 that are smaller than $1$.
    - Small fixes: start y-axis at 0, change $S$ to $E$.
    - Try decreasing record rate (maybe its the autocorrelation value)?

    __Scan past what you're supposed to read__.


# February 3
- New internal energy plot:
    - $100\times 100$
    - 100 measurements 200 sweeps apart, 2000 sweep thermalization. Wolff cluster every 5 sweeps.

![Internal Energy, revised plot](internal_energy2)
    - Mostly the same.
    - Todo:
        - Try plotting magnetic susceptibility? and compare with Table 1 [@berg1981].
        - Measure autocorrelation time (for $E$ and/or $\chi_m$). This should go up in strong coupling.

- Gradient flow:
    - after decreasing $dt$ step size, action now makes sense:
    - Implement adaptive step size if it becomes a problem.
    - Todo:
        - Begin to study other observables ($E$, $\chi_m$)

- Afterwards: work on topological susceptibility

- No meeting on Monday


![Action in flow time](flow_action)
<!--./bin/flow -b 1.0 -L 128 -o /dev/null  86.58s user 0.90s system 98% cpu 1:29.18 total-->

# February 10

- UF Acceptance!
- New computer, runs simulations faster.
 
- New result with $\chi_m$ shows discrepancy at phase transition:
    - Try removing Wolff algorithm
    - Perhaps magnetic susceptibility is related to Wolff cluster size? (see @newman1999 , p. 102) 
        - See @hasenbusch1995
    - Check out change in lattice size
    - Interesting point: the internal energy differs at a different region from the magnetic susceptibilty.

- Gradient Flow
    - Normalize $\chi_m$ in terms of $L^4$.
    - Run in strong correlation regime.


# References
