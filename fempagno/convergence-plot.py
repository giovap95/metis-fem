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
energy = 0.4999999989694232
errorNormVector1 = []
errorNormVector2 = []
errorNormVector3 = []
dofvec1 = []
dofvec2 = []
dofvec3 = []
energies = []

for h in hnumber:
    filename = 'pull-out-' + str(h)
    x, y,totdofs, energyh = inputfunction(filename)
    errorNorm1 = np.sqrt((energy-energyh)**2/energy**2)
    errorNormVector1.append(errorNorm1)
    dofvec1.append(totdofs)
    
    filename = 'pull-out-' + str(h) + '-2deg'
    x, y,totdofs, energyh = inputfunction(filename)
    errorNorm2 = np.sqrt((energy-energyh)**2/energy**2)
    errorNormVector2.append(errorNorm2)
    dofvec2.append(totdofs)
    
    filename = 'pull-out-' + str(h) + '-2deg'+'-loc-ref'
    x, y,totdofs, energyh = inputfunction(filename)
    errorNorm3 = np.sqrt((energy-energyh)**2/energy**2)
    errorNormVector3.append(errorNorm3)
    dofvec3.append(totdofs)

m1 = np.around(np.log10(errorNormVector1[4]/errorNormVector1[3])/np.log10(64/32), 2)
m2 = np.around(np.log10(errorNormVector2[4]/errorNormVector2[3])/np.log10(64/32), 2)
m3 = np.around(np.log10(errorNormVector3[4]/errorNormVector3[3])/np.log10(64/32), 2)

plt.loglog(dofvec1, errorNormVector1,'o-', dofvec2, errorNormVector2, '^-', dofvec3, errorNormVector3,'-s')   
#plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
plt.legend(['linear element, m = '+str(m1), 'quadratic element, m = '+str(m2), 'local refinement'])
plt.xlabel('Number of elements')
plt.ylabel('Energy norm')
plt.title('h-convergence')
#plt.xticks(hnumber, ['4', '8', '16', '32', '64'])

plt.show()