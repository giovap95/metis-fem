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
import time
start = time.process_time()

# import modules and specific functions
import numpy as np
import motoremesh
from boundary_conditions import BoundaryConditions
import solver

# Read mesh file from gmsh
mesh = motoremesh.GMSH('pull_out')



bcs = BoundaryConditions()
bcs.dirichlet_elements , bcs.dirichlet_nodes = bcs.find_boundary_obj(mesh,'Dirichlet')
bcs.neumann_elements , bcs.neumann_nodes = bcs.find_boundary_obj(mesh,'Neumann')
bcs.load = np.array([1000,0]).reshape((1,2)) # N/m

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


               
material_lib =           {'spring'    :             {'elastic properties' : {"Young's modulus":2.10256,
                                                                            'Poisson ratio':None},
                                                     'geometric properties': {'volumeFactor': None}},   
                  
                          'open_bar'  :             {'elastic properties' :   {"Young's modulus":2,
                                                                              'Poisson ratio':None},   
                                                    'geometric properties' : {'volumeFactor': 2}},

                          'mat1'       :            {'elastic properties' :  {"Young's modulus":100,
                                                                               'Poisson ratio' : 0.3},
                                                    'geometric properties': {'volumeFactor': 1}},

                          'matrix'  :             {'elastic properties' : {"Young's modulus":700,
                                                                            'Poisson ratio':0.3},
                                                    'geometric properties':{'volumeFactor' : 5}},
                
                
                          'fiber'  :             {'elastic properties' : {"Young's modulus":70000,
                                                                            'Poisson ratio':0.3},
                                                    'geometric properties':{'volumeFactor' : 5}},
                         }















# Solver
U,K = solver.run(mesh,bcs,material_lib,parameters)

mesh.point_data = {'Displacement':U.reshape((int(len(U)/mesh.d),mesh.d))}
meshio.write('prova2.vtk',mesh,file_format='vtk')
end = time.process_time()
print("\n...you just wasted",round(end-start,6),"seconds of your life\n \n")