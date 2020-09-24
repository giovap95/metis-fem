#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 17:51:57 2020

@author: giovanni
"""
from input import inputfunction
from plot_analytical_solutions import plot_analytical_solutions as pas
import matplotlib.pyplot as plt

x1,y1,energyh1 = inputfunction('pull-out-4-2deg-loc-ref')

x2,y2,energyh2  = inputfunction('pull-out-8-2deg-loc-ref')

x3,y3,energyh3  =  inputfunction('pull-out-16-2deg-loc-ref')

x4,y4,energyh4=  inputfunction('pull-out-32-2deg-loc-ref')

x5,y5,energyh5 = inputfunction('pull-out-64-2deg-loc-ref')

xan, yan = pas()

plt.plot(x1,y1,'-o', x2,y2,'-o',x3,y3,'-o',x4,y4,'.-',x5,y5,'r.-',xan,yan,'b--')
labels = ['4','8','16','32','64','analytical']
plt.legend(labels, title='number of elements')
plt.ylabel('Displacement')
plt.xlabel('Coordinate')
plt.title('Effect of number of 1st degree bar elements on the FE solution')
plt.show()
