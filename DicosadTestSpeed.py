from Main import main
import matplotlib.pyplot as plt

data = []

timings = range(1000,101000,1000)

complexity = []
realtime= []
simtime = []


for time in timings:
    print time
    data = main(time)
    complexity += [data[0]]
    simtime += [data[1]]
    realtime += [data[2]]

print realtime
# Three subplots sharing both x/y axes
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(timings, complexity, "b")
ax2.plot(timings, realtime, "b")
ax2.set_yscale("log")
ax3.plot(timings, simtime, "b")
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
#f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.show()