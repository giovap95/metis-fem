# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:34:55 2020

@author: giova
"""
import meshio
import sys
# specify where to look for modules
sys.path.append("modules")
sys.path.append("PRE")


# import modules and specific functions
import numpy as np
from motoremesh import Mesh
from boundary_conditions import BoundaryConditions
import solver

# Read mesh file from gmsh
gmsh = meshio.read("D:\\Documents\\GitHub\\metis-fem\\fempagno\\PRE\\composite.msh")
# Instancing classes to objects
mesh = Mesh()
element = 'triangle'

mesh.material = gmsh.get_cell_data('gmsh:physical',element)
mesh.el_def = np.ones((len(mesh.material),1)) # could be deprecated?
mesh.cds_table = gmsh.points[:,0:2] # coordinates of each node anti-clockwise

mesh.conn_table = gmsh.cells_dict[element] # nodes in each element (1 row: 1 element) anti-clockwise

mesh.elements = len(gmsh.cells_dict[element])
mesh.nodes = len(gmsh.points)
mesh.nodesperelem = 500 # deprecable
mesh.dofspernode = 2
mesh.totdofs=mesh.nodes*mesh.dofspernode


bcs = BoundaryConditions()
bcs.dirichlet_elements , bcs.dirichlet_nodes = bcs.find_boundary_obj(gmsh,'Dirichlet')
bcs.neumann_elements , bcs.neumann_nodes = bcs.find_boundary_obj(gmsh,'Neumann')
bcs.load = np.array([100,0]) 

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

                         81  :  {'element' : 'quad',
                                'elastic properties' : {"Young's modulus":70000,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 5},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'quad',
                                                        'rule':'Gauss Legendre',
                                                        'points':4}},
                         8  :  {'element'  :  'triangle',
                                'elastic properties' : {"Young's modulus":70000,
                                                        'Poisson ratio':0.3},
                                'geometric properties':{'thickness' : 5},
                                'stiffness matrix' :   {'evaluation':'numerical integration',
                                                        'domain':'triangle',
                                                        'rule':'Gauss Legendre',
                                                        'points':3}},
                         }

# Solver
U,K = solver.run(mesh,bcs,material_lib,parameters)

gmsh.point_data = {'Displacement':U.reshape((int(len(U)/2),2))}
meshio.write('prova2.vtk',gmsh,file_format='vtk')
