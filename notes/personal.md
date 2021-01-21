## Monday August 10
- Problem: After converting to Cython, observables do not indicate a phase transition at the expected value.
    - I run the a Markov chain with and without the Wolff algorithm at $m_0^2=-0.72$. The result depends on the Wolff algorithm, exhibiting the broken phase without the Wolff algorithm but the symmetric phase with the Wolff algorithm. This indicates that there is some bug in the Wolff algorithm.
    - When growing the cluster, it looks like Schaich will consider a site that has already been tested but not counted, if it is approached from a different path.
    - After implementing this change, no effect on discrepancy.
    - It looks like near the critical mass, the field is too chaotic to form meaningful clusters.
    - By looking at the list of sites in a cluster, I discovered that sometimes a point is counted twice. I think this is because neighbors are checked against the cluster BEFORE they are run. This means that a site can appear in `to_test` twice.

- Problem with Binder-Cumulant: I should be taking the ensemble average, then calculating the Binder Cumulant.
    - I implemented this change using `primary_observables` and `derived_observables`.


## Aug 25
- Issue: binder cumulant approached 2/3 very slowly below the critical mass.
- I determined that the square needs to happen after the volume-average. This fixed the issue with the binder cumulant.
- I implemented this order to compute derived quantities:
    1. calculate volume-average
    2. take any exponents
    3. take the ensemble average
    4. compute any derived quantities





## Wednesday, August 26
- Started implementing errors.
- Created "secondary observables" for powers of phi.
- Issue: BC exhibits very large errors in symmetric phase.



## Friday Aug 28

- Implemented the Jackknife method by rewriting much of the Recorder class. I am still not sure if the current implementation is the best but it seems to work for the time being. I have seperated the quantities into "primary observables" for independent quantites that are calculated for each lattice, "secondary observables" which depend on the primary observables, and "derived quantities" which depend on the ensemble averages of the secondary and primary observables.
- This method provides more realistic errors for the Binder cumulant and the susceptibility.

- The next step is a migration to the cluster. I am setting up a job script on `bora`.

## Tuesday September 1
- Had some difficulty connecting `mpi4py` on `bora`. The trick is to make sure the current `conda` implementation does not have MPI since `distutils` always check `conda` first. Also, make sure to build it from scratch where `mpicc` links to the desired MPI implementation.
- It looks like a higher number of cores does not necessarily mean faster for these small lattices.
- After running on cluster, $N=128$ looks really good.
- In order to make 12 cores work, I plan on transitioning to L=96
-Begin implementation of the gradient flow

## Thursday, September 3
- Implement the gradient flow as a new class with a single method. Use a shallow copy of the previous lattice so as not to copy the entire lattice. Reassign the `data` attribute with the evolved data.
- Implement a `hooks` parameter to the `RandomWalk.run` method, a list of functions to be run before each measurement. This is accompanied by a list of recorders that corresponds to each hook.
- ISSUE: The gradient flow makes the lattice complex. How does this affect observables?
- Start plotting the same observables as a function of flow time, not as a function of $m_0^2$.



## Tuesday, September 8
- Since the momenta in different directions are scaled equally, the inverse Fourier transform is real and should be casted as such.
- The problem in the flowed values is the ordering of the recorders. This was a simple fix which produced a much clearer plot with all observables remaining unchanged except the action, which drops sharply at $\tau=0.1$ and stays at this level.


## Friday, Septemer 25
- Begin large rewrite
    - Transition calculations to C++
    - Transition run() method to C++ so as to decrease lattice transfer of data
    - Define phi as a  vector to so that process works with nonlinear $\sigma$ model.
    

## Wednesday, October 28
- I've implemented a command line program that runs the Monte Carlo algorithm on a lattice
- It is called by a python program which should manage all the parameters, passing them as cmd line arguments to the C++ executable.
- Problems with C++ code:
    - Cluster doesn't work yet
    - MPI doesn't work yet
    - Gives wrong critical mass for $\phi^4$.

Other concerns:
    - Picking new values for $\sigma$ model, options:
        1. Using rotation matrices
        2. adding random values and then normalizing


TODO:
    - Quantitatively compare results from python and C++ codes.
        - Implement errors and use them to compare.
    -

## Friday, October 30

- Make dim and N constants of Sweeper
- Vectors are continuous! Begin shift to vectors

## Wednesday, Nov 4

- Added cluster algorithm, not successfully compiling
- Minor Issue: Progress bar stuck, doesnt move past 1% until process is done.
- Issue: results look strange at higher values of $N$
- Todo: compare python and C++ techniques
- Rewriting some of the new phi techniques

## Thursday, Jan 21
- Fix MPI Broadcasting issue
- Fixed Big Sur compatibility issue
- Implemented mpi debug system
- Generalized to sigma model
