# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 20:35:54 2020

@author: Conor
"""
from CrankNicolsonNeumann import crank_nicolson_neu
import CrankNicolsonDirichlet as cnd
import matplotlib.pyplot as plt
import numpy as np
import UsefulFunctions as uf

""" Neumann """
#def bndryA(x):
#    return 0
#
#def bndryB(x):
#    return 2
#
#def iniconds(x):
#    return 1
#
#conc = crank_nicolson_neu(0,2,3,1.5,20,200,bndryA,bndryB,iniconds)
#plt.close('all')


""" Normal diffusion dirichlet """
#x_start = 0
#x_stop = 1
#xN = 100
#tN = 100
#t_stop = 1
#dt = (t_stop)/(tN-1)
#D = 0.5
#
#def bndryA(x):
#    return 0
#def bndryB(x):
#    return 0
#def iniconds(x):
#    return np.sin(np.pi*x)
##x_start,x_stop,t_stop,D,xN,tN,bndryA,bndryB,iniconds
#conc = cnd.crank_nicolson_dir_NOPLOT(x_start,x_stop,t_stop,D,xN,tN,bndryA,bndryB,iniconds)
#
#def realsol(x,t):
#    return np.sin(np.pi*x)*np.exp(-(np.pi**2)*D*t)
#
#x = np.linspace(x_start,x_stop,xN)
#t = np.linspace(0,t_stop,tN)
#
#for i in np.arange(0,len(t),1):
#       
#    plt.plot(x,conc[:,i],'r',linewidth=3)
#    plt.plot(x[::5],realsol(x[::5],i*dt),'x',mew=5, ms=5)
#    
#    plt.ylim((0,np.amax(conc)))
#    
#    plt.pause(0.01)
#    plt.clf()
#
#plt.close('all')


x_start = 0
x_stop = 1
xN = 100
tN = 100
t_stop = 0.5
dt = (t_stop)/(tN-1)
D = 0.5
    

def bndryA(x):
    return 0
def bndryB(x):
    return 0
def iniconds(x):
    return np.sin(np.pi*x)
#x_start,x_stop,t_stop,D,xN,tN,bndryA,bndryB,iniconds
conc = cnd.crank_nicolson_dir_NOPLOT(x_start,x_stop,t_stop,D,xN,tN,bndryA,bndryB,iniconds)

def realsol(x,t):
    return np.sin(np.pi*x)*np.exp(-(np.pi**2)*D*t)

x = np.linspace(x_start,x_stop,xN)
t = np.linspace(0,t_stop,tN)

for i in np.arange(0,len(t),1):
       
    plt.plot(x,conc[:,i],'r',linewidth=3)
#    plt.plot(x[::5],realsol(x[::5],i*dt),'x',mew=5, ms=5)
    
    plt.ylim((0,np.amax(conc)))
    
    plt.pause(0.01)
    plt.clf()

plt.close('all')