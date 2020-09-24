# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:56:00 2020

@author: giova
"""
import matplotlib.pyplot as plt
import numpy as np
from input import inputfunction 

Lnumber = [20, 10, 5, 2]
error_vec = []
element_vector = []

for l in Lnumber:
    filename = 'quarter_disk' + str(l)
    error_norm, elements = inputfunction(filename)
    error_vec.append(error_norm)
    element_vector.append(elements)

m = np.log10(error_vec[3]/error_vec[2])/np.log10(element_vector[3]/element_vector[2])
plt.loglog(element_vector, error_vec,'-o')
#plt.legend('linear triangles, slope: '+str(m))
plt.xlabel('Number of elements')
plt.ylabel('Error norm')
plt.title('h-convergence')