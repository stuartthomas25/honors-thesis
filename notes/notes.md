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
-  Next week back: start gradient flow on linear phi4 model

\newpage

# References


