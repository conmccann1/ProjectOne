# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 19:26:54 2020

@author: Conor
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

x_start = 0
x_stop = 1
t_stop = 1
xN = 50
tN = 50
D = 0.25
h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)

def bndryA(x):
    return 5
def bndryB(x):
    return 1
def inicons(x):
    return 0.4

r = (D*dt)/(h**2)



x = np.linspace(0,x_stop,xN)
t = np.linspace(0,t_stop,tN)

conc = np.zeros((xN,tN))
A = np.zeros((xN-2,xN-2))
b = np.zeros((xN-2,1))

conc[0,:] = bndryA(x)
conc[-1,:] = bndryB(x)
conc[:,0] = inicons(x)
#conc[int(np.floor(xN//2)),0] = 1


A[0,0] = 2*(r+1)
A[0,1] = -r
A[-1,-1] = 2*(r+1)
A[-1,-2] = -r


for i in np.arange(1,xN-3,1):
    A[i,i] = 2*(r+1)
    A[i,i-1] = -r
    if i == xN-3:
        continue
    A[i,i+1] = -r
    

for k in np.arange(0,tN-1):
    
    b[0] = 2*conc[1,k]+r*(conc[0,k]+conc[0,k+1]+conc[2,k]-2*conc[1,k])
    
    for i in np.arange(2,xN-2):
        b[i-1] = 2*conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
        
    b[-1] = 2*conc[xN-2,k]+r*(conc[xN-1,k+1]+conc[xN-3,k]-2*conc[xN-2,k]+conc[xN-1,k])
    
    c = np.linalg.solve(A,b)
     
    for i in np.arange(0,xN-2):
        conc[i+1,k+1] = c[i]
        
def realsol(x,t):
    return np.sin(np.pi*x)*np.exp(-D*(np.pi**2)*t)

realconc = np.zeros((xN,tN))
realconc = realsol(x,t)

plt.plot(x,conc[:,1],x,realsol(x,dt))

for i in np.arange(0,len(t),1):
    plt.plot(x,conc[:,i],'b',linewidth=3)
    plt.plot(-x,conc[:,i],'b',linewidth=3)
    

#    plt.plot(x[::5],realsol(x[::5],i*dt),'x',mew=5, ms=5)
    
    plt.ylim((0,6))
    plt.pause(0.01)
    plt.clf()

plt.plot(x,conc[0,:])
#
#
#
#
