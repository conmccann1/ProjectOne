# -*- coding: utf-8 -*-
"""
Method lines for 1D diffusion equation

 -Second-order centred difference approximation for the
 spatial derivative
 
 -Euler time stepping

 -first-order in time and second-order in space.

@author: Conor
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
x_start = 0
x_stop = 1
t_stop = 2
xN = 6
tN = 12
D = 0.1
h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)

r = (D*dt)/(h**2)

#h = 0.2

initial_conds = 0.9

x = np.linspace(0,1,xN)
t = np.linspace(0,2,tN)

conc = np.zeros((len(x),len(t)))

conc[0,:] = 5
conc[-1,:] =1
conc[:,0] = initial_conds

for k in np.arange(0,len(t)-1,1):
    for i in np.arange(1,len(x)-1,1):
        conc[i,k+1] = conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
    


for i in np.arange(0,len(t),1):
    plt.plot(x,conc[:,i])
    plt.pause(0.01)
    plt.clf()
