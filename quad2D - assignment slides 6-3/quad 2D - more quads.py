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

mesh.material = [4,4,4,4]
mesh.el_def = [1,1,1,1] # could be deprecated?
mesh.cds_table = np.array([[0,0],
                           [1,0],
                           [2,0],
                           [3,0],
                           [4,0],
                           [4,1],
                           [3,1],
                           [2,1],
                           [1,1],
                           [0,1]]) # coordinates of each node anti-clockwise

mesh.conn_table = np.array([[0,1,8,9],
                            [1,2,7,8],
                            [2,3,6,7],
                            [3,4,5,6]]) # nodes in each element (1 row: 1 element) anti-clockwise
                            
mesh.elements = 4       
mesh.nodes =  10      
mesh.nodesperelem = 4       
mesh.dofspernode = 2        
mesh.totdofs=mesh.nodes*mesh.dofspernode 


bcs = BoundaryConditions()
bcs.load = np.array([[9,1],  # dof 2, forza 4
                     [11,1]]) # dof 4, forza 4

bcs.zerodisp_dof = np.array([0,1,18,19]) # dof with prescribed zero displacement

# Define parameters and the materials that will be used in the FEA

parameters = {"strain components": 3, #TODO: to be related to the spatial dimension
                                                    #strain component -> stress/strain state: 1D, 
                                                    #2D plane stress, 2D plane strain, 3D
              "solver"           : {"type": "linear"}} # try with nonlinear}
                            
material_lib =           {1  :  {'element' :  'spring',
                                'elastic properties' : {'stiffness':210},
                                'stiffness matrix':{'evaluation':'closed form'}},
    
                         2  :  {'element' :  'bar',
                                'elastic properties' :   {"Young's modulus":2},
                                'geometric properties' : {'area': 2},
                                'stiffness matrix' :     {'evaluation':'numerical integration',
                                                          'domain':'line',
                                                          'rule':'Gauss Legendre',
                                                          'points':2}},
                                
                         3  :  {'element' :  'bar',
                                'elastic properties' :  {"Young's modulus":150},
                                'geometric properties': {'area': 0.891},
                                'stiffness matrix' :    {'evaluation':'closed form'}},
                         
                         4  :  {'element' : 'quad',
                                'elastic properties' : {"Young's modulus":200,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 10},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'line',
                                                        'rule':'Gauss Legendre',
                                                        'points':2}}
                         }

# Solver
U,K = solver.run(mesh,bcs,material_lib,parameters)