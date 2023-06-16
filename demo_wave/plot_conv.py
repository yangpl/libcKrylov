import cmath
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0

#-----------------------------------------------
it, error1= np.loadtxt('bicgstab_300.txt', skiprows=0, unpack=True)
it, error2= np.loadtxt('gmres_300.txt', skiprows=0, unpack=True)


p1 = plt.plot(it, error1, label='BiCGStab')
p2 = plt.plot(it, error2, label='GMRES')
plt.xlabel('# iteration k')
plt.ylabel('$|r|$')
plt.yscale('log')
#plt.title('Normalized misfit')
plt.legend()
plt.savefig('bicgstab_gmres.png', bbox_inches='tight')
plt.show()


