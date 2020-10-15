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
from analytical_stresses import stress_analytical as str_an

def stress_recovery(mesh,U,bcs,material_lib):
    
    dimension = mesh.d # spatial dimensions
    stress_components = 3 # if dimension==2
    error_sigma = 0
    error_sigma_element = 0
    sigma_h_vec = np.zeros((mesh.elements, 3))

    
    
    for i in np.arange(mesh.elements):
        
        nodes                            =  mesh.conn_table[i]
        u                                =  np.concatenate(U[nodes,:])
        elementType                      =  mesh.elementType[i]
        evaluation, domain, rule, points = materiale.stiff_matrix_info(mesh,i)
        weights, roots                   = i_s(evaluation, domain, rule, points)
        E,ni                             = materiale.elastic_properties(mesh,material_lib,i)
        error_sigma_element = 0
        current_weight = 1
        
        D            = E / (1 - ni**2) * np.array([[1,       ni,             0],
                                                   [ni,       1,             0],
                                                   [0,        0,     .5*(1-ni)]]) # plane stress
      
        if evaluation == 'closed form':
            
            if elementType == 'triangle':
                N, B = gauss_integ.triangle(mesh, i, elementType, dimension)
                
                sigmah = (D @ (B @ u)).reshape((1, 3))
                sigma_an= str_an(mesh, i, N).reshape((3,1))
                
                
        elif evaluation == 'numerical integration':
            
            for h in range(len(weights)):
                current_roots       = roots[h]
                current_weight     = weights[h]
                dNxy, detj, N       = gauss_integ.shape_funct(mesh, i, elementType, current_roots, dimension)
                B                   = np.zeros((stress_components,len(nodes)*2))
    
                for j in range(len(dNxy[0,:])): # assemble B matrix independent of its size, automatically
                    B[0,2*j] = dNxy[0,j]
                    B[1,2*j+1] = dNxy[1,j]
                    B[2,2*j+1] = dNxy[0,j]
                    B[2,2*j] = dNxy[1,j]
                    
                sigmah = (D @ (B @ u)).reshape((1,3))
                sigma_an= str_an(mesh, i, N).reshape((3,1))
         
        sigma_h_vec[i,:] = sigmah 
        error_sigma_element += current_weight * np.linalg.norm((sigma_an - sigmah))
        error_sigma += error_sigma_element
        
    return  sigma_h_vec