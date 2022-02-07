import os
from datetime import datetime
import matplotlib
from matplotlib import pyplot as plt

os.system("cat HYDRA_{1..10}/BIRTH > BIRTHS")
os.system("stat HYDRA_{1..10}/hist.dat | grep Change | awk '{print $2,$3}' | rev | cut -c 4- | rev > DEATHS")

births = open("BIRTHS", "r").read().split("\n")
births = [i for i in births if i is not ""]
births = [datetime.strptime(i, "%a %d %b %H:%M:%S GMT %Y") for i in births]

deaths = open("DEATHS", "r").read().split("\n")
deaths = [i for i in deaths if i is not ""]
deaths = [datetime.strptime(i, "%Y-%m-%d %H:%M:%S.%f") for i in deaths]

time = [deaths[i] - births[i] for i in range(len(births))]
time = [i.total_seconds() / 3600. for i in time]

nodes = list(range(len(time)))
nodes = [i+1 for i in nodes]

speedup = [1. / (time[i] / time[0]) for i in range(len(births))]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

ax1.plot(nodes, time, color=(0.25, 0.25, 0.25))
ax1.set(xlabel = "Number of nodes (24X cores)",
        xticks = [i for i in nodes],
        ylabel = "Time [hours]",
        yticks = [2*i for i in nodes],
        title = "pop1, time reducton with multiple nodes (ASRC)")
ax1.grid(linestyle="dotted")
ax1.set_aspect(0.5, "box")

ax2.plot(nodes, speedup, color=(0.25, 0.25, 0.25), marker=".")
ax2.plot(nodes, nodes, color=(0.75, 0.75, 0.75))
ax2.set(xlabel = "Number of nodes (24X cores)",
        xticks = [i for i in nodes],
        ylabel = "Speedup",
        yticks = [i for i in nodes],
        title = "pop1, speedup with multiple nodes (ASRC)")
ax2.grid(linestyle="dotted")
ax2.set_aspect("equal", "box")

fig.tight_layout()

efficiency = [speedup[i] / nodes[i] for i in range(len(speedup))]
efficiency = [str(i)[0:4] for i in efficiency]

for x, y, eff in zip(nodes, speedup, efficiency):
    plt.annotate(eff,
            (x, y),
            color="red",
            textcoords = "offset points",
            xytext = (-5, 5),
            ha = "center")

plt.annotate("LINEAR (OPTIMAL)",
        (7,7),
        textcoords = "data",
        xytext = (6.5, 7.0),
        #ha = "center",
        rotation = 45)

fig.savefig("speedup_ASRC.png", dpi = 300)
