# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 19:41:04 2020

@author: giova
"""

import sys
import numpy as np

def integration_scheme(evaluation,domain,rule,points):
    
    if evaluation == 'numerical integration':
        
        if rule == 'Gauss Legendre':
            
            if domain == 'line':
                
                if points == 1:
                    weights = np.array([2.0])
                    roots = np.array([0.0])
                elif points == 2:
                    weights = np.array([1.0, 1.0])
                    roots = np.array([-0.5773502691896257, 
                          +0.5773502691896257])
                elif points == 3:
                    weights = np.array([0.5555555555555556, 
                           0.8888888888888888, 
                           0.5555555555555556])
                    roots = np.array([-0.7745966692414834, 
                           0.0, 
                          +0.7745966692414834])
                    
            elif domain == 'triangle':
                
                if points == 1:
                    weights = np.array([0.5])
                    roots = np.array([[0.3333333333333333,0.3333333333333333]])
                elif points == 3:
                    weights = np.array([0.16666666666666666,
                                        0.16666666666666666,
                                        0.16666666666666666])
                    roots = np.array([[0.16666666666666666,0.16666666666666666],
                                      [1.6666666666666666,0.6666666666666666],
                                      [0.6666666666666666,1.6666666666666666]])
                    
            elif domain == 'quad':
                
                if points == 1:
                    weights = np.array([4.0])
                    roots = np.array([[0.0,0.0]])
                elif points == 4:
                    weights = np.array([1.0,1.0,1.0,1.0])
                    roots = np.array([[-0.5773502691896257,+0.5773502691896257],
                                      [-0.5773502691896257,-0.5773502691896257],
                                      [+0.5773502691896257,-0.5773502691896257],
                                      [+0.5773502691896257,+0.5773502691896257]])
            
            else: 
                print('don''t know this domain')
                sys.exit()
                 
        else:
            print('dont know this rule')
            sys.exit()
    
            
    return weights,roots 

