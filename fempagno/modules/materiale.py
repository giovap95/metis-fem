# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:59:46 2020

@author: giova
"""
import sys

def elastic_properties(mesh,material_lib,i):
    
    try: 
        key = mesh.material[i]
        Young = material_lib[key]['elastic properties']["Young's modulus"]
        ni = material_lib[key]['elastic properties']['Poisson ratio']
    except KeyError:
        print('Material not defined')
        sys.exit()
        
    return Young,ni

def geometric_properties(mesh, material_lib,i):

    try: 
        key = mesh.material[i]
        volumeFactor = material_lib[key]['geometric properties']['volumeFactor']
    except KeyError:
        print('Material not defined')
            
    return volumeFactor


def stiff_matrix_info(mesh,i):

    try: 
        key = mesh.elementType[i]
        
        evaluation = mesh.element_lib[key]['stiffness matrix']['evaluation']
        domain = mesh.element_lib[key]['stiffness matrix']['domain']
        rule = mesh.element_lib[key]['stiffness matrix']['rule']
        points = mesh.element_lib[key]['stiffness matrix']['points']
        
    except KeyError:
        print('Element not defined')
        sys.exit()
    
    return evaluation, domain, rule, points

