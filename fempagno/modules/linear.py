# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:46:27 2020

@author: giova
"""
import numpy as np
from numpy.linalg import inv
import motoreFEM
import matplotlib.pyplot as plt
import sys

def linear(mesh,bcs,material_lib,parameters):
    #%% Allocating memory
    K = np.zeros(shape=(mesh.totdofs,mesh.totdofs))

    F = np.zeros(shape=(mesh.totdofs))

    # rotation matrices
    T = np.zeros((4,4))


    #%% Magic happens
print('------ Stiffness matrix assembly -------')
    for i in np.arange(mesh.elements):

        mesh.el_type(i)

        k = motoreFEM.stiffness_matrix(mesh,material_lib,parameters,T,i)

        dof = motoreFEM.locglobmap(mesh,i)

        K = motoreFEM.assembly(k,dof,K)

print('------ Applying boundary conditions ------')
    #%% Applying boundary conditions
    bcs.apply_bcs(F,K,mesh)

    #%% Solving
    #U = np.matmul(inv(K), F)

    if np.linalg.det(K)==0:
        print('\n','\n')
        print('K is singular! the geometry is not constrained. Aborting...')
        print('\n','\n')
        sys.exit()

print('------ Computing displacements ------')    
    U = inv(K)@F

    print('###########################################################')
    print('U: ',U)
    ax1 = plt.plot(mesh.cds_table[:,0],mesh.cds_table[:,1],'-o')
    Ucoord = U.reshape((mesh.nodes,2))
    m = mesh.cds_table+Ucoord
    ax2 = plt.plot(m[:,0],m[:,1],'-o')
    plt.xlim(0,2)
    plt.ylim(-1.5,2)

    plt.show()
    return U,K
