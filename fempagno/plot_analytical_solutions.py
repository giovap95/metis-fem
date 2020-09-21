#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 10:38:35 2020

@author: giovanni
"""

def plot_analytical_solutions():
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    E = 1
    A = 1
    P = 1
    #alpha1 = 0.5
    
    alpha2 = 1
    
    #alpha3 = 2
    
    #alpha4 = 4
    
    x = np.linspace(0,10,num=50)
    
    #u1 = -P/(E*A*alpha1)*np.exp(-alpha1*x)
    
    u2 = -P/(E*A*alpha2)*np.exp(-alpha2*x)
    
    #u3 = -P/(E*A*alpha3)*np.exp(-alpha3*x)
    
    #u4 = -P/(E*A*alpha4)*np.exp(-alpha4*x)

    return x,u2