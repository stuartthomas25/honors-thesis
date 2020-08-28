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
