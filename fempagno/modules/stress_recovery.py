# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:18:07 2020

@author: giova
"""
import numpy as np
from integration_scheme import integration_scheme as i_s
import materiale
import gauss_integ

def stress_recovery(mesh,U,bcs,material_lib):
    
    dimension = mesh.d # spatial dimensions
    stress_components = 3 # if dimension==2
    sigma_vec = np.zeros((mesh.elements,stress_components)) # initialize total 
    
    for i in np.arange(mesh.elements):
        
        nodes = mesh.conn_table[i]
        u = np.concatenate(U[nodes,:])
        elementType = mesh.elementType[i]
        weights, roots = i_s('numerical integration', elementType, 'Gauss Legendre', 1)
        E,ni = materiale.elastic_properties(mesh,material_lib,i)
        D = E/(1-ni**2)*np.array([[1,       ni,             0],
                                  [ni,       1,             0],
                                  [0,        0,     .5*(1-ni)]]) # plane stress
        
        for h in range(len(weights)):
            current_roots = roots[h]
            current_weights = weights[h]
            dNxy , detj = gauss_integ.shape_funct(mesh, i, elementType, current_roots, dimension)
            B = np.zeros((stress_components,len(nodes)*2))

            for j in range(len(dNxy[0,:])): # assemble B matrix independent of its size, automatically
                B[0,2*j] = dNxy[0,j]
                B[1,2*j+1] = dNxy[1,j]
                B[2,2*j+1] = dNxy[0,j]
                B[2,2*j] = dNxy[1,j]
                
            epsilon = B @ u 
            sigma = D @ epsilon 
            #sigma = sigma*detj*current_weights
        sigma_vec[i] = sigma
        Warning('Attento che funziona solo se sei in 2D con elementi lineari!!')
            
    return sigma_vec

def von_mises(sigma):
    
    sigma_eq = np.sqrt(sigma[:,0]**2+sigma[:,1]**2-sigma[:,0]*sigma[:,1]+3*sigma[:,2])
    
    return sigma_eq