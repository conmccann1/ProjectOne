# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 21:14:05 2020

@author: Conor
"""
from CrankNicolsonNeumann import crank_nicolson_neu
import CrankNicolsonDirichlet as cnd
import matplotlib.pyplot as plt
import numpy as np
import UsefulFunctions as uf

""" Set-up """
# Step sizes
x_start = 0
x_stop = 50
#xN = 300
#tN = 300
t_stop = 0.5
#dt = (t_stop)/(tN-1)
D = 1


h = 0.04291845
dt = 0.00505050505

xN = ((x_stop-x_start)/h) + 1
tN = ((t_stop)/dt) + 1

#h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)
r = (D*dt)/((h**2))

# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop,h)
t = np.arange(0,t_stop,dt)

# initialise solution matrix with given conditions
conc = np.zeros((len(x),len(t)))
#    conc[-1,:] = 2
#    conc[:,0] = 1
conc[0,:] = 0
conc[-1,:] = 0
conc[0,0] = 1/(h)

# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((len(x)-2,len(x)-2))

A[0,0] = r*(2/3)+2
A[0,1] = -r*(2/3)

A[-1,-1] = 2*(r+1)
A[-1,-2] = -r

for i in np.arange(1,len(x)-3,1):
    A[i,i] = 2*(r+1)
    A[i,i-1] = -r
    if i == xN-3:
        continue
    A[i,i+1] = -r

# set up b column vector (knowns)
b = np.zeros((len(x)-2,1))
  
""" Fill conc with approximations """
for k in np.arange(0,len(t)-1):

    b[0] = 2*conc[1,k]+r*(conc[0,k]-2*conc[1,k]+conc[2,k])

    for i in np.arange(2,len(x)-2):
        b[i-1] = 2*conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
        
    b[-1] = 2*conc[len(x)-2,k]+r*(+conc[len(x)-3,k]-2*conc[len(x)-2,k])
    b[-1] = 2*conc[len(x)-2,k]+r*(conc[len(x)-3,k]-2*conc[len(x)-2,k]+conc[len(x)-1,k]+conc[len(x)-1,k+1])

    c = np.linalg.solve(A,b)
    
    
    for i in np.arange(0,len(x)-2):
        conc[i+1,k+1] = c[i]
    
    conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]

realsol = np.zeros((len(x),len(t)))

err = np.zeros((len(x),len(t)))
errSumX = np.zeros((len(x),1))

""" Plot solutions for all time """
plt.figure()
for i in np.arange(0,len(t),1):
    err[:,i] = abs(conc[:,i]-uf.cBar(x,t[i])).T
    plt.plot(x,conc[:,i],'b',linewidth=3)
    plt.plot(x,uf.cBar(x,t[i]))
#    plt.ylim((0,cmax))
    plt.plot(np.ones(40),np.arange(0,4,0.1),'--')
    plt.plot(-1*np.ones(40),np.arange(0,4,0.1),'--')
    plt.plot(x,np.ones(len(x)),'--')
    plt.plot(-x,np.ones(len(x)),'--')
    plt.xlim((x_start,3))
    plt.pause(0.00000001)
    plt.clf()
    

errSumX = np.sum(err[:,1:],axis=1)
errSumT = np.sum(err[:,1:],axis=0)

plt.figure()
plt.subplot(211)
plt.plot(x,errSumX,linewidth=2)
plt.xlabel('x')
plt.ylabel('sum of errors')
plt.subplot(212)
plt.plot(t[1:],errSumT,linewidth=2)
plt.xlabel('t')


#plt.plot(x,conc[:,30],x,uf.cBar(x,t[30]))
    

    
