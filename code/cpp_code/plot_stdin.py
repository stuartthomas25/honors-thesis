from matplotlib import pyplot as plt
import sys

x = []
y = []
for line in sys.stdin:
    sys.stdout.write(line)
    tup = line.split(" ")
    x.append(float(tup[0]))
    y.append(float(tup[1]))

plt.plot(x,y, '.')
plt.xlabel(r"$\tau$ (flow time)")
plt.ylabel("$S$")
plt.title(r"Action in flow time, $L=128$")
plt.show()
