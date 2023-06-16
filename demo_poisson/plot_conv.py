import cmath
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0

#-----------------------------------------------
it1, error1= np.loadtxt('cg.txt', skiprows=0, unpack=True)
it2, error2= np.loadtxt('pcg.txt', skiprows=0, unpack=True)


p1 = plt.plot(it1, error1, label='CG')
p2 = plt.plot(it2, error2, label='PCG')
plt.xlabel('# iteration')
plt.ylabel('$\|r\|$')
plt.yscale('log')
#plt.title('Normalized misfit')
plt.legend()
plt.savefig('cg_pcg.png', bbox_inches='tight')
plt.show()


