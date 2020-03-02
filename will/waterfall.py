#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:27:17 2020

@author: WillMcD
"""

import matplotlib.pyplot as plt
import numpy as np


def waterfall(matrix, vector):

    fig1, ax1 = plt.subplots()
    
    wf = ax1.contourf(np.flip(np.transpose(matrix)),np.fliplr(vector))
    ax1.set_title('Waterfall')
    ax1.set_xlabel('Distance')
    ax1.set_ylabel('Time')
    
    cbar = fig1.colorbar(wf, ax=ax1)
    cbar.ax.set_ylabel('Concentration of Ca2+')