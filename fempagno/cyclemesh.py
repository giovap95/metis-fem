#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 17:51:57 2020

@author: giovanni
"""
from input import inputfunction
from plot_analytical_solutions import plot_analytical_solutions as pas
import matplotlib.pyplot as plt

x1,y1,foo = inputfunction('pull-out-4-2deg')

x2,y2,foo  = inputfunction('pull-out-8-2deg')

x3,y3,foo  =  inputfunction('pull-out-16-2deg')

x4,y4,foo  =  inputfunction('pull-out-32-2deg')

x5,y5,foo  = inputfunction('pull-out-64-2deg')

xan, yan = pas()

plt.plot(x1,y1,'-o', x2,y2,'-o',x3,y3,'-o',x4,y4,'.-',x5,y5,'r.-',xan,yan,'b--')
labels = ['4','8','16','32','64','analytical']
plt.legend(labels, title='number of elements')
plt.ylabel('Displacement')
plt.xlabel('Coordinate')
plt.title('Effect of number of 1st degree bar elements on the FE solution')
plt.show()