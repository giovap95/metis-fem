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

mesh.material = [4]
mesh.el_def = [1]
mesh.cds_table = np.array([[0,0],
                           [1,0],
                           [1,1],
                           [0,1]]) # coordinates of each node anti-clockwise

mesh.conn_table = np.array([[0,1,2,3]]) # nodes in each element (1 row: 1 element) anti-clockwise
                            
mesh.elements = 1       
mesh.nodes =  4      
mesh.nodesperelem = 4       
mesh.dofspernode = 2        
mesh.totdofs=mesh.nodes*mesh.dofspernode 

n = np.arange(mesh.nodes)
dn = mesh.dofspernode
mesh.dofs = np.array([n*dn,n*dn+1]).T.reshape((1,mesh.totdofs)) # Creating a list of the dofs

bcs = BoundaryConditions()
bcs.load = np.array([[2,4],  # dof 2, forza 4
                     [4,4]]) # dof 4, forza 4

bcs.zerodisp_dof = np.array([0,1,6,7]) # dof with prescribed zero displacement

# Define parameters and the materials that will be used in the FEA

parameters = {"strain components": 3, #TODO: to be related to the spatial dimension
                                                    #strain component -> stress/strain state: 1D, 
                                                    #2D plane stress, 2D plane strain, 3D
              "solver"           : {"type": "linear"}} # try with nonlinear}
                            
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
                                'elastic properties' :  {"Young's modulus":1.57},
                                'geometric properties': {'area': 0.891},
                                'stiffness matrix' :    {'evaluation':'closed form'}},
                         
                         4  :  {'element' : 'quad',
                                'elastic properties' : {"Young's modulus":2,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 1},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'line',
                                                        'rule':'Gauss Legendre',
                                                        'points':2}}
                         }

# Solver
U,K = solver.run(mesh,bcs,material_lib,parameters)