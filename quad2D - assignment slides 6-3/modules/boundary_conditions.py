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
        self.zerodisp_dof = None
        
    
    def apply_bcs(self,F,K,mesh):
        F[self.load[:,0]] += self.load[:,1] # prescribed nodal external forces
        
        
        F[self.zerodisp_dof] = 0 # zeroing forces on nodes with zero displacement
        
        # zeroing out zero displacement columns and rows
        K[:,self.zerodisp_dof] = 0
        K[self.zerodisp_dof,:] = 0
        K[self.zerodisp_dof,self.zerodisp_dof] = 1
