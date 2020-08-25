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





