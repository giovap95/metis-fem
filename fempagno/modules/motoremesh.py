# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:54:01 2020

@author: giova
"""

import numpy as np
import sys
import meshio


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
        
#---------------------------------------------------------------------------
# Functions below do not belong to mesh Class
#---------------------------------------------------------------------------
def el_mat(mesh,i):
    """ Returns the material of the current element, 
    as defined in the material dictionary"""
    el_mat = mesh.material[i]
    return el_mat


def el_type(mesh, i):
    """ Returns the element type of the current element""" #TODO: eliminate this function
    el_type = mesh.elementType[i]
    
    if el_type!=0 and el_type!=1:
                
        print('\n','Element', i, 'ERROR! Element type not recognised')
        sys.exit()
                
    return el_type

def coordinates(mesh,i):
    rows = mesh.conn_table[i]
    cds = mesh.points[rows]
              
    return cds    

    
def NodesInElement(mesh,i):
    
    NodesInElement=mesh.conn_table[i]
    
    return NodesInElement
    
    
def get_key(my_dict,val): 
    """ Function to return key for any value. """
    
    # This function returns the key if the first item in the array value 
    # of a dictionary is equal to val. If my_dict contains 
    # 'Fixed': array([667,   0]), get_key(my_dict,667) returns Fixed
    
    for key, value in my_dict.items(): 
         if val == value[0]: 
             return key 
    
    print("\n value",val,"doesn't exist as \'key\': array([value, 0]) in\n", my_dict)
    sys.exit()

      
        
def GMSH(mesh_file):
    
    sys.path.append("PRE")
    # create a mesh object
    mesh = meshio.read("D:/Documents/GitHub/metis-fem/fempagno/PRE/"+mesh_file+".msh")    
    # check if the mesh object contains attributes needed by pyFEM
    # - pyFEM_MeshAttributes is a list of all the mesh attributes needed by pyFEM
    # - we are going to reuse the attribute points and add the other attribute from pyFEM_MeshAttributes
    pyFEM_MeshAttributes = ["d", "dofsNode", "elements", "elementMaterialTag", "elementType", "points"]

    for attribute in pyFEM_MeshAttributes:
        if attribute in dir(mesh):
            if attribute == "points":
                pass
            else:
                print("Error: meshio already contains the attribute",attribute)
                print("       ...do something!")
                sys.exit()

    # add the missing attributes from pyFEM_MeshAttributes
    
    # Note: it is assumed that the mesh is two-dimensional and that the
    # domain is dicretized with triangular elements and that there are
    # two degrees of freedom per node (i.e., this is a plain equilibrium problem)
    
    mesh.elements = 0
    mesh.nodes = len(mesh.points)
    mesh.dofspernode = 2
    mesh.totdofs=mesh.nodes*mesh.dofspernode
    mesh.d = 2 
    mesh.dofsNode = 2 
    mesh.conn_table = []
    mesh.material = []
    mesh.el_def = []
    mesh.elementType = []
    mesh.material = []
    meshing = False
    
    
    quad = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['quad']
        quad = True
    except KeyError:
        # print("No quadrilateral elements in mesh")
        pass
    
    triangle = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['triangle']
        triangle = True
    except KeyError:
        # print("No triangular elements in mesh")
        pass

    if quad:
        meshing = True
        quads = len(mesh.cell_data_dict["gmsh:physical"]["quad"])
        mesh.elements += quads
        for t in range(quads):
            mesh.conn_table.append(mesh.cells_dict["quad"][t])
            materialTag=mesh.cell_data_dict["gmsh:physical"]["quad"][t]
            # we assume that a physical surface in 2D is only used to identify 
            # elements with the same material property.
            # GMSH identifies a physical group by a tag and a name. 
            # Tags are stores in cell_data_dict for each element.
            # Tags and names are linked in field_data              
            # The function get_key returns the name (=key) for a given tag            
            key = get_key(mesh.field_data, materialTag)
            mesh.material.append(key)            
            mesh.elementType.append('quad')
        
    if triangle: 
        meshing = True
        triangles = len(mesh.cell_data_dict["gmsh:physical"]["triangle"])
        mesh.elements += triangles
        for t in range(triangles):
            mesh.conn_table.append(mesh.cells_dict["triangle"][t])
            materialTag=mesh.cell_data_dict["gmsh:physical"]["triangle"][t]
            # we assume that a physical surface in 2D is only used to identify 
            # elements with the same material property.
            # GMSH identifies a physical group by a tag and a name. 
            # Tags are stores in cell_data_dict for each element.
            # Tags and names are linked in field_data              
            # The function get_key returns the name (=key) for a given tag            
            key = get_key(mesh.field_data, materialTag)
            mesh.material.append(key)            
            mesh.elementType.append('triangle')
    
    if not meshing:
        print("something went wrong: could not extract mesh data")
        sys.exit()
    
    mesh.points = mesh.points[:, 0:mesh.d]     #resize to the number of spatial dimensions in the problem
    # TODO: ...check that all the necessary attributes have been defined in a correct manner


# library of the possible elements
    mesh.element_lib =  { 'spring'    :    {'stiffness matrix'  :  {'evaluation' : 'closed form',
                                                                    'domain'     : None,
                                                                    'rule'       : None,
                                                                    'points'     : None}},
                                                             
                          'bar'      :    {'stiffness matrix'  :   {'evaluation' : 'numerical integration',
                                                                    'domain'     : 'line',
                                                                    'rule'       : 'Gauss Legendre',
                                                                    'points'     :  2}},
                          
                          'triangle'  :    {'stiffness matrix'  :  {'evaluation' : 'numerical integration',
                                                                    'domain'     : 'triangle',
                                                                    'rule'       : 'Gauss Legendre',
                                                                    'points'     : 1}},
                                                
                          'quad'      :    {'stiffness matrix'  : {'evaluation'  : 'numerical integration',
                                                                   'domain'      : 'quad',
                                                                   'rule'        : 'Gauss Legendre',
                                                                   'points'      : 4}}
                          }



    return mesh
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

