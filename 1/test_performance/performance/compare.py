import matplotlib.pyplot as plt
import numpy as np

data0 = np.genfromtxt('loop_noflag.csv', delimiter=',')
data1 = np.genfromtxt('loop_flag.csv', delimiter=',')

#print(data)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
axes[0].plot(data0[:,1],data0[:,0], '-o', ms=2)
axes[1].plot(data1[:,1],data1[:,0], '-o', ms=2)
axes[0].set_ylim((0,.0055))
axes[1].set_ylim((0,.0055))
fig.subplots_adjust(wspace=0.2)
axes[1].set_yticklabels([])
axes[0].set_title("Loop no flag")
axes[1].set_title("Loop -O3 flag")

#plt.show()
plt.savefig('cputimes.png')
