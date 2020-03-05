# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:44:58 2020

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
x_start = -5
x_stop = 5
#xN = 300
#tN = 300
t_stop = 0.5
#dt = (t_stop)/(tN-1)
D = 17


h = 0.01
dt = 0.001

xN = ((x_stop-x_start)/h) + 1
tN = ((t_stop)/dt) + 1

#h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)
r = (D*dt)/((h**2))

# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop+h,h)
t = np.arange(0,t_stop,dt)

siteloc = np.arange(1/h,2*x_stop/h,1/h)


conc = np.zeros((len(x),len(t)))
#    conc[-1,:] = 2
#    conc[:,0] = 1
conc[0,:] = 0

# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((len(x),len(x)))

conc[-1,:] = 0
conc[int(len(x)/2),0] = 430

A[0,0] = 1 + (2)*r
A[0,1] = -(2)*r

A[-1,-1] = 1 + (2)*r
A[-1,-2] = -(2)*r

for i in np.arange(0,len(x),1):
    A[i,i] = 2*(r)+1
    A[i,i-1] = -r
    if i == len(x)-1:
        break
    A[i,i+1] = -r
    
# set up b column vector (knowns)
b = np.zeros((len(x),1))


siteL = 1
siteR = 1
sparktimes = np.zeros(x_stop*2)

for k in np.arange(0,len(t)-1,1):
        
    for i in np.arange(0,len(x)):
        b[i] = conc[i,k]
    c = np.linalg.solve(A,b)
    print(k)
    for i in np.arange(0,len(x)):
        conc[i,k+1] = c[i]
    
    conc[0,k+1] = c[1]
    conc[-1,k+1] = c[-2]

#std = np.diff(sparktimes)
#sm1 = conc[400,:]
#s1 = conc[600,:]
#
#sdif = abs(sm1-s1)

err = np.zeros((len(x),len(t)))
errSumX = np.zeros((len(x),1))

plt.figure()
for i in np.arange(0,len(t),1):
    err[:,i] = abs(conc[:,i]-uf.cBar(x,t[i])).T
    plt.plot(x,conc[:,i],'b',linewidth=3)
#    plt.plot(x[::10],uf.cBar(x[::10],t[i]),'ro',mew=3,ms=3)
###    plt.ylim((0,cmax))
#    
#    plt.plot(np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(-1*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(2*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(-2*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(3*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(-3*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(4*np.ones(40),np.arange(0,4,0.1),'--')
#    plt.plot(-4*np.ones(40),np.arange(0,4,0.1),'--')
#    
#    plt.plot(x,np.ones(len(x)),'--')
#    plt.plot(-x,np.ones(len(x)),'--')
#    plt.xlim((-5,5))
    plt.pause(0.00000001)
#    
    if i == len(t)- 1:
        break
    
    plt.clf()
##    
#    
errSumX = np.sum(err[:,1:],axis=1)
errSumT = np.sum(err[:,1:],axis=0)

plt.figure()
plt.subplot(211)
plt.title('1-Norm Error for Implicit Method')
plt.plot(x,errSumX,linewidth=2)
plt.xlabel('x')
plt.ylabel('sum of errors')
plt.subplot(212)
plt.plot(t[1:],errSumT,linewidth=2)
plt.ylabel('sum of errors')
plt.xlabel('t')
#plt.close('all')