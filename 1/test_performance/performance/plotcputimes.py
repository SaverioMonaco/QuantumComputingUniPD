import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('cputimes.csv', delimiter=',')

#print(data)

plt.plot(data[:,1],data[:,0], '-o', ms=3)
#plt.show()
plt.savefig('cputimes.png')
