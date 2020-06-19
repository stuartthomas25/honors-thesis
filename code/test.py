from timeit import timeit
from matplotlib import pyplot as plt
import math
import numpy as np

times = list(range(1,50))
reassign = []
no_reassign = []
def numpy():
    np.exp(-1000.0)

def math():
    math.exp(-1000.0)

for t in times:
    reassign.append( timeit(math, globals=dict( t=t)))
    no_reassign.append( timeit(numpy, globals=dict(t=t)))

plt.plot(times,reassign, label="math")
plt.plot(times,no_reassign, label="numpy")
plt.legend()
plt.show()
