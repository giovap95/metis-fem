# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 12:44:55 2020

@author: giova
"""

import numpy as np
import sys

def bar(mesh,young,area,weights,roots,cds,parameters,i):
    """ Compute the stiffness matrix through Gauss_Legendre integration of a 1D bar
    of the first order in a 1D domain"""

    length = np.sqrt((cds[1][0]-cds[0][0])**2+(cds[1][1]-cds[0][1])**2) # finds the length through pythagora's theorem
    
    jac = length/2
    B = np.array([-1/length,1/length]).reshape(1,2)
    
    k = jac * np.sum(weights) * B.T@B
    k = k*young*area/length
    return k


def triangle(mesh,E,ni,t,weights,roots,cdsreal,parameters):
    """ Compute the stiffness matrix of a generic 2D triangular element 
    in a 2D-ONLY domain using the isoparametric approach"""
    
    k = np.zeros((6,6))
    
    D = E/(1-ni**2)*np.array([[1,       ni,             0], # TODO: implement D in the parameters input
                              [ni,       1,             0],
                              [0,        0,     .5*(1-ni)]]) # plane stress
    dN = np.array([[1,0,-1],
                   [0,1,-1]]).T
    
    x = cdsreal[:,0]
    y = cdsreal[:,1]
    jac = np.zeros((2,2))
    jac[0,0] = np.sum(dN[:,0]*x)
    jac[0,1] = np.sum(dN[:,0]*y)
    jac[1,0] = np.sum(dN[:,1]*x)
    jac[1,1] = np.sum(dN[:,1]*y) #TODO: validate
    
    det = np.linalg.det(jac)
    
    dNxy = np.zeros((3,2))
    dNxy = dN@np.linalg.inv(jac)
    
    B = np.zeros((3,6))
    B = np.array([[dNxy[0,0],0,dNxy[0,1],0,dNxy[0,2],0],
                  [0,dNxy[1,0],0,dNxy[1,1],0,dNxy[1,2]],
                  [dNxy[1,0],dNxy[0,0],dNxy[1,1],dNxy[0,1],dNxy[1,2],dNxy[0,2]]])
    
    k = B.T@D@B*det*t
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
            w = weights[j]*weights[h]
            dN = np.zeros((4,2))

            
            # For a parametric quadrilateral element the derivative shape function matrix is:
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
            print('csi,eta old \n',csi,eta)
            print('dNxy old \n',dNxy)
            print('jac old \n',jac)
    return k


def shape_funct(mesh, i, el_type, roots, dim):
    nodesinelement = len(mesh.NodesInElement(i))
    N = np.zeros((1,nodesinelement))
    dN = np.zeros((dim,nodesinelement))
    jac = np.zeros((dim,dim))
    el_coord = mesh.coordinates(i)   # form: [x1,y1],
                                     #       [x2,y2],... CRUCIAL
    
    #choose the correct parametric shape function
    if el_type == 'triangle':
        s = roots[0]
        r = roots[1]
        N = np.array([1-r-s , s , r])
        
        dN = np.array([[-1,1,0],
                       [-1,0,1]])
        
    elif el_type == 'quad':
        csi = roots[0]
        eta = roots[1]
        N = np.array([.25*(1-csi)*(1-eta) , .25*(1+csi)*(1-eta) , .25*(1+csi)*(1+eta) , .25*(1-csi)*(1+eta)])

        dN   =  np.array([[-.25*(1-eta) ,  .25*(1-eta)  , .25*(1+eta) , -.25*(1+eta)],
                          [-.25*(1-csi) , -.25*(1+csi) ,  .25*(1+csi) ,  .25*(1-csi)]])
    else:
        print('WARNING: No procedure coded for this element')
        sys.exit()   
        
        
    jac = dN @ el_coord
    dNxy = np.linalg.inv(jac) @ dN
    detj = np.linalg.det(jac)
    return dNxy, detj
    
    
                          

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    