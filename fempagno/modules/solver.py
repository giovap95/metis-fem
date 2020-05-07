# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:56:24 2020

@author: giova
"""

def run(mesh, bcs, material_lib, parameters):

    solverType = parameters["solver"]["type"]
    
    exec("from "+solverType+" import "+solverType)
    
    U,K = eval (solverType+"(mesh, bcs, material_lib, parameters)")    
    
    # ...you might want to have a safer approach to eval:
    # https://www.geeksforgeeks.org/eval-in-python/
    
    return U,K