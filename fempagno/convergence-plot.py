# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:45:56 2020

@author: giova
"""

from input import inputfunction 
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# Cycle through element size h
hnumber = [4, 8, 16, 32, 64]
hsize = [10/4, 10/8, 10/16,  10/32, 10/64]
energy = 0.5
errorNormVector1 = []
errorNormVector2 = []
energies = []

for h in hnumber:
    filename = 'pull-out-' + str(h)
    x, y, energyh = inputfunction(filename)
    energies.append(energyh)
    errorNorm = np.sqrt((energy-energyh)/energy)
    errorNormVector1.append(errorNorm)
    
    filename = 'pull-out-' + str(h) + '-2deg'
    x, y, energyh = inputfunction(filename)
    errorNorm = np.sqrt((energy-energyh)/energy)
    errorNormVector2.append(errorNorm)


plt.loglog(hnumber, errorNormVector1,'o-', hnumber, errorNormVector2, '^-')   
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
plt.legend(['linear element', 'quadratic element'])
plt.xlabel('Number of elements')
plt.ylabel('Energy norm')
plt.title('h-convergence')
#plt.xticks(hnumber, ['4', '8', '16', '32', '64'])

plt.show()