# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:15:40 2020

@author: giova
"""
import numpy as np
import motoreFEM
import scipy.sparse as sps


def nonlinear(mesh,bcs,material_lib,parameters):
    
    u_bar = 0.5 
    T = 10
    Dt = 0.01 
    n = 3 
    toll = 1e-6
    E = 1e6
    A = 100 
    L = .5 
    
    F = np.zeros(shape=(mesh.totdofs))
    u = np.zeros((3,1))
    fi = np.zeros((3,1))
    fe = np.zeros((3,1))
    h = np.zeros((4, 1))
    
    c = []
    r = []
    data = []
    t = 0 
    
    while t<T:
        t = t+Dt
        duu = np.array([0, u_bar*Dt])
        j = 0 
        q = np.inf
        while q>toll and j<10:
            j = j+1;
            # assembly of tangent stiffness matrix
            for i in range(mesh.elements):
                k = motoreFEM.stiffness_matrix(mesh,material_lib,parameters,T,i)
        
                dof = motoreFEM.locglobmap(mesh,i)
        
                r_new, c_new, data_new = motoreFEM.assembly(k,dof)
                
                r.extend(r_new)
                c.extend(c_new)
                data.extend(data_new)
                
            K = sps.coo_matrix((data,(r,c)),shape=(mesh.totdofs,mesh.totdofs))
            #K = K.tocsr()
            K = K.toarray()
            ### END OF ASSEMBLY ###
            #%% Applying boundary conditions
            
            iiu = bcs.find_dofs(mesh,mesh.cell_sets_dict['Dirichlet']['vertex']) # constrained dofs
            iif = np.setdiff1d(bcs.find_dofs(mesh,mesh.conn_table), iiu) # unconstrained dofs
            Kfu = K[iif,iiu]
            Kff = K[iif,iif]
            fif = fi[iif]
            fef = fe[iif]
            
            ### SOLVER ###
            duf = Kff/(fef-fif-Kfu@duu)
            u[iiu] = u[iiu] + duu
            u[iif] = u[iif] + duf
            duu = np.zeros((2,1))
            
            ### INTERNAL STRESS###
            #fi = #int(B.t*sigma)dOmega
            
            ### COMPUTE REACTION FORCES ###
            fe[iiu] = fi[iiu]
            
            ### RESIDUAL ###
            r = np.norm(fi-fe)/(1+np.norm(fe))
            
            
            
