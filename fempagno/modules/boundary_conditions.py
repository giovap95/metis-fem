# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:33:51 2020

@author: giova
"""
import numpy as np

class BoundaryConditions:
    
    def __init__(self):
        self.load = None
        self.disp = None
        self.dirichlet_nodes = None
        self.neumann_nodes = None
        
    def find_dofs(self,mesh,nodes):
        dn = mesh.dofspernode
        dofs = np.array([nodes * dn,
                          nodes * dn + 1])
        dofs = np.concatenate(dofs.T)
        return dofs
    
    
    def apply_bcs(self,F,K,mesh):
        
        # Find dofs where Neumann conditions are enforced
        neumann_dofs = self.find_dofs(mesh,self.neumann_nodes)
        # Apply forces on neumann dofs
        F[neumann_dofs] += self.load # prescribed nodal external forces
        
        
        # Find dofs where Dirichlet conditions are enforced
        dirichlet_dofs = self.find_dofs(mesh,self.dirichlet_nodes)
        F[dirichlet_dofs] = 0 # zeroing forces on nodes with zero displacement
        # zeroing out zero displacement columns and rows
        K[:,dirichlet_dofs] = 0
        K[dirichlet_dofs,:] = 0
        K[dirichlet_dofs,dirichlet_dofs] = 1
