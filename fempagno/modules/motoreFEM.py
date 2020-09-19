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
import motoremesh
import scipy.sparse as sps
## Global stiffness matrix
#%% generate k_i
def stiffness_matrix(mesh,material_lib,parameters,T,i):
    evaluation, domain, rule, points = materiale.stiff_matrix_info(mesh,i)

    elMat = mesh.material[i]
    elType = mesh.elementType[i]

    d = parameters["strain components"] #TODO: move out of the iterating stiffness_matrix function
    dim = parameters["spatial dimensions"]
#%%
    if evaluation == 'closed form':

        if elType == 'spring': # check if element is a spring
            spring_stiffness = materiale.elastic_properties(mesh,material_lib,i)
            k = spring_stiffness * np.array([[1, -1],
                                             [-1, 1]])


        elif elType == 'bar': #check if element is a bar
            young = materiale.elastic_properties(mesh,material_lib,i)
            area = materiale.geometric_properties(mesh,material_lib,i)
            cds = motoremesh.coordinates(mesh,i)
            length = np.abs(np.diff(cds)) # finds the length
            k = (young * area)/length * np.array([[1,-1],
                                                  [-1, 1]])

        else:
            print("don't know what to do")
            sys.exit()



#%%
    if evaluation == 'numerical integration':
        
        k = 0
        E,ni = materiale.elastic_properties(mesh,material_lib,i)
        t = materiale.geometric_properties(mesh,material_lib,i)
        weights,roots = i_s(evaluation,domain,rule,points)
        D = E/(1-ni**2)*np.array([[1,       ni,             0],
                                  [ni,       1,             0],
                                  [0,        0,     .5*(1-ni)]]) # plane stress

        for h in range(len(weights)):
            current_roots = roots[h]
            dNxy , detj = gauss_integ.shape_funct(mesh, i, elType, current_roots, dim)
            B = np.zeros((d,len(motoremesh.NodesInElement(mesh,i))*2))

            for j in range(len(dNxy[0,:])): # assemble B matrix independent of its size, automatically
                B[0,2*j] = dNxy[0,j]
                B[1,2*j+1] = dNxy[1,j]
                B[2,2*j+1] = dNxy[0,j]
                B[2,2*j] = dNxy[1,j]
            k += B.T@D@B*t*detj*weights[h]

    return k


#%% global mapping
def locglobmap(mesh,i): # extract dof of current element

    n = np.array(mesh.conn_table[i])
    dn = mesh.dofspernode
    dof = np.array([n*dn]) #np.zeros(mesh.nodesperelem, dtype=np.uint8)
    dof = np.concatenate(dof.T)
    return dof

#%% assembly of global stiffness matrix
def assembly(k_mat,cur_dof):
   
   prod = [(x,y) for x in cur_dof for y in cur_dof]
   r = [x for (x,y) in prod]
   c = [y for (x,y) in prod]
   data = np.concatenate(k_mat)
   #K[np.ix_(cur_dof,cur_dof)] += k_mat


   return r, c, data
