# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:59:46 2020

@author: giova
"""
import sys

def elastic_properties(elMat,material_lib,i):
    
    key = elMat

    if material_lib[key]['element'] == 'spring':
        springStiffness =material_lib[key]['elastic properties']['stiffness']            
        return springStiffness
    
    elif material_lib[key]['element'] == 'bar':            
        E = material_lib[key]['elastic properties']["Young's modulus"]            
        return E
    
    elif material_lib[key]['element'] == 'quad':
        E = material_lib[key]['elastic properties']["Young's modulus"]
        ni = material_lib[key]['elastic properties']['Poisson ratio']
        return E,ni
    else:
        print('Material not defined')
        sys.exit()


def geometric_properties(elMat,material_lib,i):
    
    key = elMat
    
    if material_lib[key]['element'] == 'bar':            
        A = material_lib[key]['geometric properties']['area']            
        return A
    
    elif material_lib[key]['element'] == 'quad':
        t = material_lib[key]['geometric properties']['thickness']
        return t
    else:
        print('geometry of element',i,'not recognised')
        sys.exit()


def stiff_matrix_info(mesh,material_lib,i):
    
    key = mesh.material[i]
    
    element = material_lib[key]['element']
    print('-------------------------------------------------------------')
    print('element',i,'is a ',element)
    evaluation = material_lib[key]['stiffness matrix']['evaluation']
    print('evaluation: ',evaluation)
    
    if evaluation == 'numerical integration':
        domain = material_lib[key]['stiffness matrix']['domain']
        print('domain: ',domain)
        
        rule = material_lib[key]['stiffness matrix']['rule']
        print('rule: ',rule)
        
        points = material_lib[key]['stiffness matrix']['points']
        print('points: ',points)
        return element, evaluation, domain, rule, points
    elif evaluation == 'closed form':
        return element, evaluation, None, None, None
        

    
    