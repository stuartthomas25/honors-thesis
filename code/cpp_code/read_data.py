import numpy as np
from matplotlib import pyplot as plt
import csv

actions = []
phibars = []
m02s = [-0.80, -0.76, -0.72, -0.68, -0.64];

with open("actions.csv", 'r') as datafile:
    csv_data = csv.reader(datafile)
    for lines in csv_data:
        actions.append(float(lines[0]))
        phibars.append(float(lines[1]))

plt.plot(m02s, phibars, '.', label=r"$\bar{\phi}$")
plt.legend()
plt.show()



