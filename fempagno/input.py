# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:34:55 2020

@author: giova
"""

import sys
# specify where to look for modules
sys.path.append("modules")


# import modules and specific functions
import numpy as np
from motoremesh import Mesh
from boundary_conditions import BoundaryConditions
import solver


# Instancing classes to objects
mesh = Mesh()

mesh.material = [5]
mesh.el_def = [1] # could be deprecated?
mesh.cds_table = np.array([[0,0],
                           [100,0],
                           [100,100]]) # coordinates of each node anti-clockwise

mesh.conn_table = np.array([[0,1,2]]) # nodes in each element (1 row: 1 element) anti-clockwise
                            
mesh.elements = 1      
mesh.nodes = 3      
mesh.nodesperelem = 500 # deprecable      
mesh.dofspernode = 2        
mesh.totdofs=mesh.nodes*mesh.dofspernode 


bcs = BoundaryConditions()
bcs.load = np.array([[5,1000]]) # dof x, forza y (in newton N)

bcs.zerodisp_dof = np.array([0,1,2]) # dof with prescribed zero displacement

# Define parameters and the materials that will be used in the FEA

parameters = {"strain components": 3, #TODO: to be related to the spatial dimension
                                                    #strain component -> stress/strain state: 1D, 
                                                    #2D plane stress, 2D plane strain, 3D
              "solver"           : {"type": "linear"},
              "spatial dimensions": 2} # try with nonlinear}
 

# Unit of measure:
# E = MPa
# F = N
# A = mm^2
# t = mm
                           
material_lib =           {1  :  {'element' :  'spring',
                                'elastic properties' : {'stiffness':2.10256},
                                'stiffness matrix':{'evaluation':'closed form'}},
    
                         2  :  {'element' :  'bar',
                                'elastic properties' :   {"Young's modulus":2},
                                'geometric properties' : {'area': 2},
                                'stiffness matrix' :     {'evaluation':'numerical integration',
                                                          'domain':'line',
                                                          'rule':'Gauss Legendre',
                                                          'points':2}},
                                
                         3  :  {'element' :  'bar',
                                'elastic properties' :  {"Young's modulus":100},
                                'geometric properties': {'area': 1},
                                'stiffness matrix' :    {'evaluation':'closed form'}},
                         
                         4  :  {'element' : 'quad',
                                'elastic properties' : {"Young's modulus":70000,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 5},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'quad',
                                                        'rule':'Gauss Legendre',
                                                        'points':4}},
                         5  :  {'element'  :  'triangle',
                                'elastic properties' : {"Young's modulus":70000000,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 5},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'triangle',
                                                        'rule':'Gauss Legendre',
                                                        'points':3}}, #TODO: change points' domain to 4 to avoid "for" loops (also for quads)
                         }

# Solver
U,K = solver.run(mesh,bcs,material_lib,parameters)
