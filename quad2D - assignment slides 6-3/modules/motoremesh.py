# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:54:01 2020

@author: giova
"""

import numpy as np
import sys

# Creates a mesh class 
class Mesh:
    
    def __init__(self):
        
        self.el_def = None
        self.material = None
        self.conn_table = None
        self.cds_table = None
        self.elements = None # globdof.shape[0]        
        self.nodes = None # max(max(globdof[:,-1]),max(globdof[:,-2]))+1       
        self.nodesperelem = None         
        self.dofspernode = None    
        self.totdofs = None
        self.d = None # spatial dimensions
        
    def coordinates(self,i):
        row = self.conn_table[i]
        cds = self.cds_table[row]
                  
        return cds    
    
        
    def NodesInElement(self,numelem):
        
        NodesInElement=self.conn_table[numelem]
        
        return NodesInElement
    
    
    def el_type(self, i):
        """ Returns the element type of the current element"""
        el_type = self.el_def[i]
        
        if el_type!=0 and el_type!=1:
                    
            print('\n','Element', i, 'ERROR! Element type not recognised')
            sys.exit()
                    
        return el_type
        
    def el_mat(self,i):
        """ Returns the material of the current element, 
        as defined in the material dictionary"""
        el_mat = self.material[i]
        return el_mat
        
        