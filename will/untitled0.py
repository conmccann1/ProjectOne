#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:21:59 2020

@author: WillMcD
"""

import numpy as np
import matplotlib as plt
import mpmath as mp
import math

t = np.linspace(0.001,5,1000)
y = [mp.nsum(lambda n:math.exp(-n/4*x)/math.sqrt(4*np.pi*x*n),[1,mp.inf]) for x in t]

plt.figure()
plt.plot(t,y)