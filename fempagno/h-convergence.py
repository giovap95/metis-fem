# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:56:00 2020

@author: giova
"""
import matplotlib.pyplot as plt
import numpy as np
from input import inputfunction 

Lnumber = [50, 20, 10]
triangle_error_vec = []
quad_error_vec = []
triangle_element_vector = []
quad_element_vector = []

for l in Lnumber:
    filename = 'quarter_disk' + str(l)
    error_norm, elements = inputfunction(filename)
    triangle_error_vec.append(error_norm)
    triangle_element_vector.append(elements)
    
    filename = 'quarter_disk' + '-q-'+str(l)
    error_norm, elements = inputfunction(filename)
    quad_error_vec.append(error_norm)
    quad_element_vector.append(elements)
    
m = np.log10(triangle_error_vec[2]/triangle_error_vec[1])/np.log10(triangle_element_vector[2]/triangle_element_vector[1])
plt.loglog(triangle_element_vector, triangle_error_vec,'-o', quad_element_vector, quad_error_vec,'^-')
plt.legend(['triangles','quads'])
plt.xlabel('Number of elements')
plt.ylabel('Error norm')
plt.title('h-convergence')