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
    import matplotlib.pyplot as plt
    
    # import modules and specific functions
    import numpy as np
    import motoremesh
    from boundary_conditions import BoundaryConditions
    import solver
    
    # Read mesh file from gmsh
    mesh = motoremesh.GMSH(filename)
    
    
    
    bcs = BoundaryConditions()
    bcs.dirichlet_elements , bcs.dirichlet_nodes = bcs.find_boundary_obj(mesh,'end')
    bcs.neumann_elements , bcs.neumann_nodes = bcs.find_boundary_obj(mesh,'start')
    bcs.load = 1 # N/mm
    
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
    
    
    
    material_lib =           {'spring'    :             {'elastic properties' : {"Young's modulus":2000,
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
    
    
                              'AISI 316'  :             {'elastic properties' : {"Young's modulus":1,
                                                                                'Poisson ratio':0.3,
                                                                                "matrix stiffness":2},
                                                        'geometric properties':{'volumeFactor' : 1}},
                             }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Solver
    U,K = solver.run(mesh,bcs,material_lib,parameters)
    energyh = .5 * U.T @ K @ U
    U = np.concatenate((U[[0]], U[2:], U[[1]]))
    
    
    # plot
    x = mesh.points[:,0]
    x.sort()
    y = U
    #plt.plot(x,y, '.-')
    
    # Writing data
    
    #point_data = {'Displacement':U}
    #meshio.write_points_cells('prova2.vtk', mesh.points, mesh.cells, point_data=point_data, cell_data = mesh.cell_data)
    #meshio.write('prova2.vtk', mesh, file_format='vtk', cell_data=cell_data)
    
    end = time.process_time()
    print("\n...you just wasted",round(end-start,6),"seconds of your life\n \n")
    return x,y,energyh