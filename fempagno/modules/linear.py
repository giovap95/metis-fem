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
import time
from tqdm import trange
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve


def linear(mesh,bcs,material_lib,parameters):
    #%% Allocating memory
    # K = np.zeros(shape=(mesh.totdofs,mesh.totdofs))
    
    r = []
    c = []
    data = [] # storing variables for assembly

    F = np.zeros(shape=(mesh.totdofs))

    # rotation matrices
    T = np.zeros((4,4))


    #%% Magic happens

    print('------ Stiffness matrix assembly -------')

    for i in trange(mesh.elements):
        
        #mesh.el_type(i)

        k = motoreFEM.stiffness_matrix(mesh,material_lib,parameters,T,i)
        
        dof = motoreFEM.locglobmap(mesh,i)

        r_new, c_new, data_new = motoreFEM.assembly(k,dof)
        
        r.extend(r_new)
        c.extend(c_new)
        data.extend(data_new)
        
    K = sps.coo_matrix((data,(r,c)),shape=(mesh.totdofs,mesh.totdofs))
    K = K.tocsr()

    print('------ Applying boundary conditions ------')
    #%% Applying boundary conditions
    bcs.apply_bcs(F,K,mesh)

    #%% Solving
    #U = np.matmul(inv(K), F)

    #if np.linalg.det(K)==0:
    #    print('\n','\n')
    #    print('K is singular! the geometry is not constrained. Aborting...')
     #   print('\n','\n')
     #   sys.exit()


    print('------ Computing displacements ------')

    #U = inv(K)@F
    U = spsolve(K, F)

    print('###########################################################')
    print('first ten U entries: ',U[0:10])

    return U,K
