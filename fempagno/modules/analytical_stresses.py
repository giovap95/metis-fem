# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 14:39:45 2020

@author: giova
"""
import motoremesh
import numpy as np

def stress_analytical(mesh,i,N):
### CALCOLO SOLUZIONE ANALITICA ###
# Ricavo le coordinate dell'elemento corrente
    coordinates = motoremesh.coordinates(mesh,i)
    x = coordinates[:,0]
    y = coordinates[:,1]
    x_root = N @ x
    y_root = N @ y 
    a = 25 #mm, radius of the hole
    P = 1 #N/mm2 load at infinite
    r = np.sqrt(x_root**2 + y_root**2) # distance from the centre
    theta = np.arccos(x_root/r) # angle in radians
    
    sigma_an_xx = P * (1 - a**2/r**2*(3/2*np.cos(2*theta)+np.cos(4*theta))+(3*a**4)/(2*r**4)*np.cos(4*theta))
    sigma_an_yy = P * (-a**2/r**2*(1/2*np.cos(2*theta)-np.cos(4*theta))-(3*a**4)/(2*r**4)*np.cos(4*theta))
    sigma_an_xy = P * (-a**2/r**2*(1/2*np.sin(2*theta)+np.sin(4*theta))+(3*a**4)/(2*r**4)*np.sin(4*theta))
    
    sigma_an = np.array([[sigma_an_xx, sigma_an_yy, sigma_an_xy]]).reshape((3, 1))
    
    return sigma_an