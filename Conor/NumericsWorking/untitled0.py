# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:22:48 2020

@author: Conor
"""

from CrankNicolsonNeumann import crank_nicolson_neu
import CrankNicolsonDirichlet as cnd
import matplotlib.pyplot as plt
import numpy as np
import UsefulFunctions as uf

a = (1/((2*np.pi*np.e)**0.5))

""" Set-up """
# Step sizes
x_start = 0
x_stop = 50
#xN = 300
#tN = 300
t_stop = 1
#dt = (t_stop)/(tN-1)
D = 1


h = 0.01291845
dt = 0.00505050505

xN = ((x_stop-x_start)/h) + 1
tN = ((t_stop)/dt) + 1

#h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)
r = (D*dt)/((h**2))

# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop,h)
t = np.arange(0,t_stop,dt)

conc = np.zeros((len(x),len(t)))
#    conc[-1,:] = 2
#    conc[:,0] = 1
conc[0,:] = 0

# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((len(x)-2,len(x)-2))

conc[-1,:] = 0
conc[500,0] = 32
#A[0,0] = r*(2/3)+2
#A[0,1] = -r*(2/3)
A[0,0] = 1 + (2/3)*r
A[0,1] = -(2/3)*r

A[-1,-1] = 2*(r)+1
A[-1,-2] = -r

for i in np.arange(1,len(x)-3,1):
    A[i,i] = 2*(r)+1
    A[i,i-1] = -r
    if i == xN-3:
        continue
    A[i,i+1] = -r
    
# set up b column vector (knowns)
b = np.zeros((len(x)-2,1))
  
for k in np.arange(0,len(t)-1,1):
    
    for i in np.arange(1,len(x)-1):
        b[i-1] = conc[i,k]
    c = np.linalg.solve(A,b)
    
   
    
    for i in np.arange(0,len(x)-2):
        conc[i+1,k+1] = c[i] 
    
    conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]


plt.figure()
for i in np.arange(0,len(t),1):
#    err[:,i] = abs(conc[:,i]-uf.cBar(x,t[i])).T
    plt.plot(x,conc[:,i],'b',linewidth=3)
#    plt.plot(x[::3],uf.cBar(x[::3],t[i]),'rx',mew=5,ms=5)
##    plt.ylim((0,cmax))
#    plt.plot(np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(-1*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(x,np.ones(len(x)),'--')
#    plt.plot(-x,np.ones(len(x)),'--')
#    plt.xlim((x_start,3))
    plt.pause(0.00000001)
    plt.clf()
#    