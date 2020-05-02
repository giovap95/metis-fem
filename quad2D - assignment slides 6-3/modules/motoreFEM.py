# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 17:54:24 2020

@author: giova
"""
import sys
# specify where to look for modules
sys.path.append("fempagno")
import numpy as np
from integration_scheme import integration_scheme as i_s
import gauss_integ
import materiale
## Global stiffness matrix

# generate k_i
def stiffness_matrix(mesh,material_lib,parameters,T,i):
    element, evaluation, domain, rule, points = materiale.stiff_matrix_info(mesh,material_lib,i)
    
    elMat = mesh.el_mat(i)
    
    if evaluation == 'closed form':
        
        if element == 'spring': # check if element is a spring
            spring_stiffness = materiale.elastic_properties(elMat,material_lib,i)
            k = spring_stiffness * np.array([[1, -1],
                                             [-1, 1]])
            
        
        elif element == 'bar': #check if element is a bar
            young = materiale.elastic_properties(elMat,material_lib,i)
            area = materiale.geometric_properties(elMat,material_lib,i)
            cds = mesh.coordinates(i)
            length = np.sqrt((cds[1][0]-cds[0][0])**2+(cds[1][1]-cds[0][1])**2) # finds the length through pythagora's theorem
            print('length of element',i,' = ',length)
            k = (young * area)/length * np.array([[1,0, -1,0],
                                                    [0,0,0,0],
                                                  [-1,0, 1,0],
                                                  [0,0,0,0]])
            # Rotation to match the real axes
            c = (cds[1][0]-cds[0][0])/length
            s = (cds[1][1]-cds[0][1])/length
            T = np.array([[c,s,0,0],
                          [-s,c,0,0],
                          [0,0,c,s],
                          [0,0,-s,c]])
            k = T.T@k@T
            
            
            
        else:
            print("don't know what to do")
            sys.exit()
    
    if evaluation == 'numerical integration':
        
        if element == 'bar':
            ## NB does'nt work if bar is not parallel to x-axis
            young = materiale.elastic_properties(elMat,material_lib,i)
            area = materiale.geometric_properties(elMat,material_lib,i)
            cds = mesh.coordinates(i)                             # ATTENTION! BETTER DEFINITION OF CDS IS NEEDED
            weights,roots = i_s(evaluation,domain,rule,points)
            
            k = gauss_integ.bar(mesh,young,area,weights,roots,cds,parameters,i)
            
        elif element == 'quad': 
            E,ni = materiale.elastic_properties(elMat,material_lib,i)
            t = materiale.geometric_properties(elMat,material_lib,i)
            weights,roots = i_s(evaluation,domain,rule,points)
            cdsreal = mesh.coordinates(i)
            
            k = gauss_integ.quad(mesh,E,ni,t,weights,roots,cdsreal,parameters)
    
        else:
            print("don't know what to do")
            sys.exit()
            
    return k

# global mapping
def locglobmap(mesh,i): # extract dof of current element
    
    n = mesh.conn_table[i]
    dn = mesh.dofspernode
    dof = np.array([n*dn,n*dn+1]) #np.zeros(mesh.nodesperelem, dtype=np.uint8)
    dof = np.concatenate(dof.T)
    return dof


# assembly of global stiffness matrix
def assembly(k_mat,cur_dof,K):

   K[np.ix_(cur_dof,cur_dof)] += k_mat   
  
   
   return K
