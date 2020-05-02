# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 12:44:55 2020

@author: giova
"""

import numpy as np

def bar(mesh,young,area,weights,roots,cds,parameters,i):
    """ Compute the stiffness matrix through Gauss_Legendre integration of a 1D bar
    of the first order in a 1D domain"""

    length = np.sqrt((cds[1][0]-cds[0][0])**2+(cds[1][1]-cds[0][1])**2) # finds the length through pythagora's theorem
    
    jac = length/2
    B = np.array([-1/length,1/length]).reshape(1,2)
    
    k = jac * np.sum(weights) * B.T@B
    k = k*young*area/length
    return k




def quad(mesh,E,ni,t,weights,roots,cdsreal,parameters):
    """ Compute the stiffness matrix of a generic 2D quadrilateral element in a 
    2D-ONLY domain using the isoparametric approach"""
    
    cdsnat = np.array([[-1,-1],[1,-1],[1,1],[-1,1]])
    
    k = np.zeros((8,8))
    D = E/(1-ni**2)*np.array([[1,       ni,             0],
                              [ni,       1,             0],
                              [0,        0,     .5*(1-ni)]]) # plane stress
    
    # two nested loops are needed because the integration in a planar domain is done
    # over roots.size*2 points (4 points in case of a 2D element)
    for h in range(2):
        
        for j in range(2):
            csi = roots[h]
            eta = roots[j]
            w = weights[j]
            dN = np.zeros((4,2))
            ni=0.3
            
            # For a parametric quadrilateral element the shape function matrix is:
            for i in range(4):
               dN[i,0]=.25*(cdsnat[i,0])*(1+cdsnat[i,1]*eta)
               dN[i,1]=.25*(1+csi*cdsnat[i,0])*(cdsnat[i,1])
               
            
            # Compute the jacobian matrix
            jac = np.zeros((2,2))
            
            jac[0,0] = np.sum(dN[:,0]*cdsreal[:,0])
            jac[0,1] = np.sum(dN[:,0]*cdsreal[:,1])
            jac[1,0] = np.sum(dN[:,1]*cdsreal[:,0])
            jac[1,1] = np.sum(dN[:,1]*cdsreal[:,1])
            
            # determinant of the jacobian matrix
            det = np.linalg.det(jac)
            
            dNxy = np.zeros((4,2))
            
            # The shape function's derivative is mapped from (csi,eta) to (x,y)
            dNxy = dN@np.linalg.inv(jac)
            
            B = np.zeros((3,8))
            
            # Matrix of shape functions derivative
            B = np.array([[dNxy[0,0],0,dNxy[1,0],0,dNxy[2,0],0,dNxy[3,0],0],
                          [0,dNxy[0,1],0,dNxy[1,1],0,dNxy[2,1],0,dNxy[3,1]],
                          [dNxy[0,1],dNxy[0,0],dNxy[1,1],dNxy[1,0],dNxy[2,1],dNxy[2,0],dNxy[3,1],dNxy[3,0]]])
            
            
            # Stiffness matrix
            k = k + B.T@D@B*det*w*t
    return k