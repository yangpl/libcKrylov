import cmath
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0

#-----------------------------------------------
it1, error1= np.loadtxt('cgnr.txt', skiprows=0, unpack=True)
it2, error2= np.loadtxt('cgne.txt', skiprows=0, unpack=True)


p1 = plt.plot(it1, error1, label='CGNR')
p2 = plt.plot(it2, error2, label='CGNE')
plt.xlabel('# iteration')
plt.ylabel('$\|r\|$')
plt.yscale('log')
#plt.title('Normalized misfit')
plt.legend()
plt.savefig('cgnr_cgne.png', bbox_inches='tight')
plt.show()


