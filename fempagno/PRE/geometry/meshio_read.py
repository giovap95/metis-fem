#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 11:13:42 2020

GMSH and meshio

@author: angelo.simone@unipd.it
"""

# meshio: https://github.com/nschloe/meshio

# install the package meshio from the Spyder console using
# conda install -c conda-forge meshio 

# This script has been tested with latest (28 March 2020) GMSH version: 4.5.5

# Do yourself a favour and use GMSH dark interface: tools, options, general, use dark interface


# You could use the GMSH api and obtain all the information you need directly from GMSH. 
# Instructions are on
# https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py
# with a quick start on
# https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api_start/
# ...or use meshio. 
# You could also modify the geo file on the fly and generate a new mesh with GMSH
# (this could be useful in the mesh refinement study in the second assignment)

# First thing to do is to get familiar with GMSH. A very good introduction 
# can be found on https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/

# Focus on:    
# 1) https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_firstmesh/
# 2) https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_advance/
# 3) https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_physical_vs_elementary/    
    
import meshio
import numpy as np
import os
from numpy.lib.arraysetops import unique

#...for pretty printing
np.set_printoptions(precision=3)

# -----------------------------------------------------------------------------

def pause():
    programPause = input("Press <ENTER> to continue.")

# -----------------------------------------------------------------------------

print("\n meshio version", meshio.__version__,"\n")

print("\n\n\n\n----------- PART 1 -----------------\n\n\n\n\n")

# Here you could have instructions that write/modify the geo file. You could
# start from an existing geo file or write one using the GMSH API
# (see https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api_start/)
# On-the-fly modifications to the geo file can be handy if you perform a 
# mesh refinement study and would like to check the effect of changes 
# you made to your FEM code.

# GMSH can be executed from python by using a system call as done below.

# Give the name of the geo file without extension 
# ...in this way it will be easier to assign the extension to the output file
# square1 comes from https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_physical_vs_elementary/
mesh_file = "square1"

# run GMSH agains the input .geo file and generate a .msh file
os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")
# this instruction is equivalent to 
# os.system("gmsh -2 square1.geo -o square1.msh")
# Check sections 3.2 and 3.3 of the GMSH manual (https://gmsh.info/doc/texinfo/gmsh.html)
# for command-line options

# ...the "correct" instructions should be the following:

# try:
#     os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")
#     gmsh = True
# except OSError:
#     print("Gmsh needs to be installed for this example.")
#     gmsh = False
    
# and then execute the rest only if gmsh is true:

# if gmsh:
#     do something
# else:
#     exit


# import the discretization as an object
mesh = meshio.read(mesh_file+".msh")

# ...what's inside the object?
# to print all the attributes (data and methods) of an object use the dir() function
#print(dir(mesh))

# ...you actually need those without underscores
attributes = [attr for attr in dir(mesh) 
              if not attr.startswith('__')]
print("\n ---------->  list of attributes\n", attributes)
 
# print attributes
for p in attributes:
    print("\n \n----------> attribute: ", p)
    string="mesh."+p
    print("what is it?",type(eval(string)))
    print(eval(string))

print("\n \n ...or you can check directly by \"printing\" the object\n")
print(mesh)
    
# to visualize the data, follow the instructions on 
# check https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_physical_vs_elementary/

pause()

print("\n\n\n\n----------- PART 2 -----------------\n\n\n\n\n")

# To identify groups of elements, GMSH uses physical tags. These tags can be used
# to identify nodes on the boundary (to enforce boundary conditions) or region 
# in a domain (to assign material properties):
# https://gmsh.info/doc/texinfo/gmsh.html#Elementary-entities-vs-physical-groups


# The geo file square2 is taken from section 3 of 
# https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/basics_physical_vs_elementary/
# and introduces a physical line.

# The geo file square3 uses a friendlier way to identify boundaries by means of labels

# The geo file bimarial adds an inclusion to the domain

# The code below uses bimaterial.geo to extract relevant data from the 
# output produced by GMSH. 

# Before you use these instructions in your code, take a look at the meshes 
# produced by GMSH with the geo files that come with this script

# -----------------------------------------------------------------------------

mesh_file = "bimaterial"
os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")
mesh = meshio.read(mesh_file+".msh")

# -----------------------------------------------------------------------------

attributes = [attr for attr in dir(mesh) 
              if not attr.startswith('__')]
for p in attributes:
    print("\n \n----------> attribute: ", p)
    string="mesh."+p
    print("what is it?",type(eval(string)))
    print(eval(string))
   
# -----------------------------------------------------------------------------

print("\n \n----------> the object contains the following data: ")
print(mesh)

# -----------------------------------------------------------------------------

print("\n \n----------> extracted attributes: points (coordinates [x y z]) \n")

# Since GMSH counts also the "central" node used to generate the arc, you'll find
# an extra node in the GSMH output and statistics. 
# The actual mesh nodes, those attached to elements, are in the coordinate array.

points = mesh.points
print(points)
print("\n size of coordinate array:",points.shape)
print("\n the number of points in the mesh can be extracted from the number of rows using len(points)")
print(" number of mesh points:",len(points))

# -----------------------------------------------------------------------------

print("\n \n----------> extracted attributes: cells (connectivity tables) \n")

# To read the mesh generated with recombined.geo, you need to add a check on
# cell.type == "quad"

# For simplicity, use only one type of domain element in your mesh (i.e., quads or triangle)
# Line elements are used (in 2D) to indentify the boundary

triangle_cells = []
line_cells = []
for cell in mesh.cells:
    if cell.type == "triangle":
        if len(triangle_cells) == 0:
            triangle_cells = cell.data
        else:
            triangle_cells = np.vstack([triangle_cells, cell.data])
    elif  cell.type == "line":
        if len(line_cells) == 0:
            line_cells = cell.data
        else:
            line_cells = np.vstack([line_cells, cell.data])

print("\n line_cells\n",line_cells)
print("\n number of line_cells:",len(line_cells))
print("\n triangle_cells\n",triangle_cells)
print("\n number of triangle_cells:",len(triangle_cells))

print("\n ...you could also use mesh.cells_dict that has the advantage of")
print("having all the connectivities grouped per entity type (vertex, line,")
print("triangle...) in a dictionary.")
line_cells1 = mesh.cells_dict["line"]
triangle_cells1 = mesh.cells_dict["triangle"]

print("\n line_cells1\n",line_cells1)
print("\n number of line_cells1:",len(line_cells1))
print("\n triangle_cells1\n",triangle_cells1)
print("\n number of triangle_cells1:",len(triangle_cells1))

print("\n or you can use a built-in method to get the same information")
triangle_cells2 = mesh.get_cells_type("triangle")
print("\n triangle_cells2 (with method) \n",triangle_cells2)
print("\n number of triangle_cells2:",len(triangle_cells2))

# -----------------------------------------------------------------------------

print("\n \n----------> extracted attributes: cell_data_dict (element tags) \n")

print("\n this lists the physical tag of each entity")

vertex_data = []
line_data = []
triangle_data = []
for key in mesh.cell_data_dict["gmsh:physical"].keys():
    if key == "vertex":
        vertex_data = mesh.cell_data_dict["gmsh:physical"][key]
    elif key == "line":
        line_data = mesh.cell_data_dict["gmsh:physical"][key]
    elif key == "triangle":
        triangle_data = mesh.cell_data_dict["gmsh:physical"][key]

print("\n cell data for vertex, size",len(vertex_data), "\n",vertex_data)
print("\n cell data for line, size",len(line_data), "\n",line_data)
print("\n cell data for triangle, size",len(triangle_data), "\n",triangle_data)

print("\n or extract info by accessing cell_data_dict as a dictionary")
triangle_data1 = mesh.cell_data_dict["gmsh:physical"]["triangle"]
print("\n cell data1 for triangle, size",len(triangle_data1), "\n",triangle_data1)

print("\n ...or use a method")
print("\n cell data for triangle\n",mesh.get_cell_data("gmsh:physical","triangle"))

# -----------------------------------------------------------------------------

print("\n \n----------> extracted attributes: field_data (physical regions) \n")

# field data are stored in a dictionary and link element tags to physical regions
#
# mesh.field_data for bimaterial.geo:
# {'bottom left corner': array([5, 0]), 'Dirichlet': array([3, 1]), 'Neumann': array([4, 1]), 'matrix': array([1, 2]), 'inclusion': array([2, 2])}
# Each key corresponds to the physical group name in the .geo file. 
# Each value is expressed as [<physical_tag>, <dimension>] (a point has dimension 0)
# In the composite mesh, 'bottom left corner': array([5, 0]) means that 
# 'bottom left corner' is a point (dimension = 0) and has tag = 5

print("mesh.field_data:\n", mesh.field_data)

print("\n this can be used to retrieve the tag of a physical region from its label")
print("\n tag for label \"Dirichlet\":", mesh.field_data["Dirichlet"][0])
print("\n tag for label \"matrix\":", mesh.field_data["matrix"][0])

# -----------------------------------------------------------------------------

print("\n \n----------> extracted attributes: mesh.cell_sets_dict (entity number by physical region label)")
print("\n Note that each entity type has its own list: vertex, line, triangle...\n ")
print("mesh.cell_sets_dict",mesh.cell_sets_dict)
print("\n points in set with label Points\n",mesh.cell_sets_dict["Points"])
print("\n elements in matrix\n",mesh.cell_sets_dict["matrix"])
print("\n line elements with Dirichlet BCs\n",mesh.cell_sets_dict["Dirichlet"]["line"])
print("\n line elements with Neumann BCs\n",mesh.cell_sets_dict["Neumann"]["line"])

# -----------------------------------------------------------------------------

print("\n \n----------> ...in practice: mesh data \n")

#total number of elements of type triangle
triangles = mesh.cells_dict["triangle"].shape[0]

for e in range(triangles):
    print("\n\nTriangle #",e,"with material tag", triangle_data1[e])
    print("\n coordinates\n", points[triangle_cells1[e]])

# -----------------------------------------------------------------------------

print("\n \n----------> ...in practice: boundary nodes (from boundary elements) \n")

print("\n line elements with Dirichlet tag\n",mesh.cell_sets_dict["Dirichlet"]["line"])

name = "Dirichlet"
tag = mesh.field_data[name][0]
dim = mesh.field_data[name][1]

if dim == 0:
    # array containing indices of elements in the boundary
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["vertex"] == tag)[0]
    # array containing indices of nodes in the boundary
    nodes = unique(mesh.cells_dict["vertex"][on_boundary])
elif dim == 1:
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
    nodes = unique(mesh.cells_dict["line"][on_boundary])

print("\n nodes related to tag Dirichlet\n", nodes)
for n in nodes:
    print("\nnode #",n,"@", points[n])
print("\n\n")


print("\n node entries in dictionary with tag \"Points\"",mesh.cell_sets_dict["Points"]["vertex"])

name = "Points"
tag = mesh.field_data[name][0]
dim = mesh.field_data[name][1]

if dim == 0:
    # array containing indices of elements in the boundary
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["vertex"] == tag)[0]
    # array containing indices of nodes in the boundary
    nodes = unique(mesh.cells_dict["vertex"][on_boundary])
elif dim == 1:
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
    nodes = unique(mesh.cells_dict["line"][on_boundary])

print("\n ...which correspond to the mesh node numbers (zero-based):", nodes)

for n in nodes:
    print("\nnode #",n,"@", points[n])
print("\n\n")

# -----------------------------------------------------------------------------

pause()

print("\n\n\n\n----------- PART 3 -----------------\n\n\n\n\n")

print("meshio can also write a vtk file")

meshio.write(mesh_file+".vtk", mesh)

# you can find examples on how to pass specific information such as the displacement field on
# https://www.litianyi.me/2019/10/27/paraview-meshio-plugin/ (see under meshio)
