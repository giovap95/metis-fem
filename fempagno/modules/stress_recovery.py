# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:18:07 2020

@author: giova
"""
import numpy as np
from integration_scheme import integration_scheme as i_s
import materiale
import gauss_integ
import motoremesh

def stress_recovery(mesh,U,bcs,material_lib):
    
    dimension = mesh.d # spatial dimensions
    stress_components = 3 # if dimension==2
    sigmah_vec = np.zeros((mesh.elements,stress_components)) # initialize total 
    sigma_an_vec = np.zeros((mesh.elements, stress_components))
    error_sigma = 0
    
    for i in np.arange(mesh.elements):
        
        nodes = mesh.conn_table[i]
        u = np.concatenate(U[nodes,:])
        elementType = mesh.elementType[i]
        weights, roots = i_s('numerical integration', elementType, 'Gauss Legendre', 1)
        E,ni = materiale.elastic_properties(mesh,material_lib,i)
        D = E/((1+ni)*(1-2*ni))*np.array([[1-ni,       ni,             0],
                                          [ni,       1-ni,             0],
                                          [0,        0,     .5*(1-2*ni)]]) # plane strain
        
        for h in range(len(weights)):
            current_roots = roots[h]
            current_weights = weights[h]
            dNxy , detj, N = gauss_integ.shape_funct(mesh, i, elementType, current_roots, dimension)
            B = np.zeros((stress_components,len(nodes)*2))

            for j in range(len(dNxy[0,:])): # assemble B matrix independent of its size, automatically
                B[0,2*j] = dNxy[0,j]
                B[1,2*j+1] = dNxy[1,j]
                B[2,2*j+1] = dNxy[0,j]
                B[2,2*j] = dNxy[1,j]
                
            epsilon = B @ u * detj 
            sigmah = D @ epsilon 
            
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
            sigmah = sigmah.reshape((3, 1))
            
            error_sigma_element = (sigma_an - sigmah).T @ (sigma_an - sigmah)
            
        sigmah_vec[i] = sigmah.reshape((1,3))
        sigma_an_vec[i] = np.array([sigma_an_xx, sigma_an_yy, sigma_an_xy])
        
        error_sigma += error_sigma_element
    print('Attento che lo stress recovery funziona solo se sei in 2D con elementi lineari!!')
            
    ### CALCOLO SOLUZIONE ANALITICA ###
    # Ricavo le coordinate dell'elemento corrente
    return sigmah_vec, sigma_an_vec, error_sigma

def von_mises(sigma):
    
    sigma_eq = np.sqrt(sigma[:,0]**2+sigma[:,1]**2-sigma[:,0]*sigma[:,1]+3*sigma[:,2]**2)
    
    return sigma_eq