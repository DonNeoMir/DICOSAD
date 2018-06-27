from Main import main
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

n = 100

components = []
histos = []
for i in range(n):
    data = main(0)
    components += [data]
    polymers = map(len,data)
    x,y = np.histogram(polymers,bins=max(polymers), range=(1,max(polymers)))
    histos += [x]

print histos

a = np.zeros((n,len(max(histos, key=len))))
for i in range(len(histos)):
    for j in range(len(histos[i])):
        a[i][j] = histos[i][j]
    
print a

plt.bar(range(1,len(max(histos, key=len))+1), np.mean(a,axis=0), yerr = stats.sem(a))

plt.show()