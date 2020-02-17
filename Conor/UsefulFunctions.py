# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 20:25:07 2020

Useful Functions

---Contains--- 
 cBar(x,t)

@author: Conor
"""
import numpy as np

def cBar(x,t):
    a = (1/((2*np.pi*np.e)**0.5))
    return np.exp(-(x**2) / (4*t)) / ( ((4*t*np.pi*(a**2)))**0.5 )

#h = 0.04291845
#dt = 0.00505050505
