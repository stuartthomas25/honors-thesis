from matplotlib import pyplot as plt
import sys

x = []
y = []
y2 = []
for line in sys.stdin:
    sys.stdout.write(line)
    tup = line.split(" ")
    x.append(float(tup[0]))
    y.append(float(tup[1]))
    y2.append(float(tup[2]))
plt.figure()
plt.plot(x,y, '.')
plt.xlabel(r"$\tau$ (flow time)")
plt.ylabel("$S$")
plt.title(r"Action in flow time, $L=128$")
plt.figure()
plt.plot(x,y2, '.')
plt.xlabel(r"$\tau$ (flow time)")
plt.ylabel("$\chi_m$")
plt.title(r"Magnetic susceptibility in flow time")
plt.show()
