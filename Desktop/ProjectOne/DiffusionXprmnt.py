# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:56:05 2020

@author: Conor
"""
import matplotlib.pyplot as plt
import numpy as np
from CrankyDiffusion1D import crankydiffusion1D
from CrankyDifNeu import crankydiffusionNEU
from CrankyDiffusion1D import crankydiffusionNOPLOT
#from CrankyNeuTEST import crankydiffusionNEUtest
#from CrankDifNeuTest2 import crankdiffusionNEUtest1


""" Normal Diffusion Dirichlet BCs """
D = 0.25
x_start = 0
x_stop = 1
xN = 50
tN = 50
t_stop = 1
dt = (t_stop)/(tN-1)

def bndryA(x):
    return 0
def bndryB(x):
    return 0
def iniconds(x):
    return np.sin(np.pi*x)

conc = crankydiffusionNOPLOT(x_start,x_stop,t_stop,bndryA,bndryB,iniconds,D,xN,tN)

def realsol(x,t):
    return np.sin(np.pi*x)*np.exp(-(np.pi**2)*D*t)

x = np.linspace(x_start,x_stop,xN)
t = np.linspace(0,t_stop,tN)

for i in np.arange(0,len(t),1):
       
    plt.plot(x,conc[:,i],'r',linewidth=3)
    plt.plot(x[::5],realsol(x[::5],i*dt),'x',mew=5, ms=5)
    
    plt.ylim((0,np.amax(conc)))
    
    plt.pause(0.01)
    plt.clf()

plt.close('all')






""" INSULATED BAR """
###### x_start,x_stop,t_stop,D,xN,tN,spike
a = (1/((2*np.pi*np.e)**0.5)) 
conc = crankydiffusionNEU(0,2,3,1.5,20,100,1/a)
plt.close('all')



""" CLOSEST TO PAPER BUT SPIKE OVER RANGE OF X AT T = 0 """
a = (1/((2*np.pi*np.e)**0.5)) 
conc = crankydiffusionNEUtest(0,5,0.5,1,50,100,1/a)
plt.close('all')



""" MOVING DIRICHLET """
def bndryA(x):
    return a-x/38
def bndryB(x):
    return 0
def iniconds(x):
    return 0

crankydiffusion1D(0,5,1,bndryA,bndryB,iniconds,0.5,50,50)
plt.close('all')


""" NOTES """
# neumann at either end???

#
#a = (1/((2*np.pi*np.e)**0.5)) 
#crankdiffusionNEUtest1(0,5,3,1.5,20,100,1/a)
#plt.close('all')