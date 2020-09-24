# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:34:55 2020

@author: giova
"""

def inputfunction(filename):
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
    from stress_recovery import stress_recovery, von_mises
    
    # Read mesh file from gmsh
    mesh = motoremesh.GMSH(filename)
    
    bcs = BoundaryConditions()
    bcs.dirichlet_elementsLeft , bcs.dirichlet_nodesLeft = bcs.find_boundary_obj(mesh,'DirichletLeft')
    bcs.dirichlet_elementsBottom , bcs.dirichlet_nodesBottom = bcs.find_boundary_obj(mesh,'DirichletBottom')
    bcs.neumann_elements , bcs.neumann_nodes = bcs.find_boundary_obj(mesh,'Neumann')
    bcs.load = np.array([1,0]).reshape((1,2)) # N/mm
    
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
    
    
                              'AISI 316'  :             {'elastic properties' : {"Young's modulus":210000,
                                                                                'Poisson ratio':0.29},
                                                        'geometric properties':{'volumeFactor' : 1}},
                             }
    
    
    # Solver
    U,K = solver.run(mesh,bcs,material_lib,parameters)
    
    # Postprocessing
    U = U.reshape((int(len(U)/mesh.d),mesh.d)) # reshaping U vector to match spatial dimensions (u_x, u_y, u_z)
    sigmah, sigma, sigma_error = stress_recovery(mesh,U,bcs,material_lib)
    
    error_norm = np.linalg.norm(sigma-sigmah)/np.linalg.norm(sigma)
    
    
    sigma_vm = von_mises(sigmah)
    
    
    
    
    # Writing data
    cells = {}
    try:
        cells['triangle'] = mesh.cells_dict['triangle']
    except KeyError:
        pass
    
    try:
        cells['quad'] = mesh.cells_dict['quad']
    except KeyError:
        pass
    
    point_data = {'Displacement':U}
    cell_data = {'Stress':sigmah,
                 'Von-Mises':sigma_vm}
    meshio.write_points_cells('prova2.vtk', mesh.points, cells, point_data = point_data, cell_data = cell_data)    
    
    
    end = time.process_time()
    print("\n...you just wasted",round(end-start,6),"seconds of your life\n \n")
    return error_norm, mesh.elements