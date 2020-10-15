# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:29:01 2020

@author: giova
"""

import matplotlib.pyplot as plt
from input import inputfunction 

Lnumber = [10, 5, 3, 2]
error_vec = []
elements_vec = []
triangle = True
quad = True


for l in Lnumber:
    if triangle:

        filename = 'quarter_disk-with fillet' + str(l)
        error_sigma, elements = inputfunction(filename)
        elements_vec.append(elements)
        plt.plot(elements, error_sigma,'ob')        


        
    if quad:
        
        filename = 'quarter_disk-with filletq' + str(l)
        error_sigma, elements = inputfunction(filename)
        elements_vec.append(elements)
        plt.plot(elements, error_sigma,'or')         
        

plt.title('convergenza norma dello stress')
plt.xlabel('numero di elementi')
plt.ylabel('norma dello stress')

plt.legend('triangoli', 'quadrilateri')
plt.show()
