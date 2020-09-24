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
alphavec = [0.5, 1,  2, 4]
energy = 0.4999999989694232
errorNormVector1 = []
errorNormVector2 = []
energies = []
label= []

i = 0
erroralpha = np.zeros((4, 5))

for alpha in alphavec:
    
    for h in hnumber:
        #filename = 'pull-out-' + str(h)
        #x, y, energyh = inputfunction(filename, alpha)
        #errorNorm1 = np.sqrt((energy-energyh)**2/energy**2)
        #errorNormVector1.append(errorNorm1)
        
        filename = 'pull-out-' + str(h) + '-2deg'
        x, y, energyh = inputfunction(filename,alpha)
        errorNorm2 = np.sqrt((energy-energyh)**2/energy**2)
        errorNormVector2.append(errorNorm2)

    
    erroralpha[i,:] = errorNormVector2
    plt.loglog(hnumber,erroralpha[i,:])
    label.append(str(alpha))
    i = i+1
    errorNormVector2 = []
#m1 = np.around(np.log10(errorNormVector1[4]/errorNormVector1[3])/np.log10(64/32), 2)
#m2 = np.around(np.log10(errorNormVector2[4]/errorNormVector2[3])/np.log10(64/32), 2)

    

#plt.loglog(hnumber, errorNormVector1,'o-', hnumber, errorNormVector2, '^-')   
#plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
plt.legend(label)
plt.xlabel('Number of elements')
plt.ylabel('Energy norm')
plt.title('h-convergence')
#plt.xticks(hnumber, ['4', '8', '16', '32', '64'])

plt.show()