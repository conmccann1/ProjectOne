# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:35:22 2020

@author: Conor
"""
import matplotlib.pyplot as plt
import numpy as np


def crankydiffusion1D(x_start,x_stop,t_stop,bndryA,bndryB,iniconds,D,xN,tN):
        
    """ Set-up """
    # Step sizes
    h = (x_stop-x_start)/(xN-1)
    dt = (t_stop)/(tN-1) 
    r = (D*dt)/(h**2)
    
    # space and time vecotors, used to fill bndrys and iniconds
    x = np.linspace(x_start,x_stop,xN)
    t = np.linspace(0,t_stop,tN)
    
    # initialise solution matrix with given conditions
    conc = np.zeros((xN,tN))
    conc[0,:] = bndryA(x)
    conc[-1,:] = bndryB(x)
    conc[:,0] = iniconds(x)

    # initialise 'A' matrix (unknowns) to be solved at each time step 
    A = np.zeros((xN-2,xN-2))
    
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
    
    # set up b column vector (knowns)
    b = np.zeros((xN-2,1))
    
    
    """ Fill conc with approximations """
    for k in np.arange(0,tN-1):
    
        b[0] = 2*conc[1,k]+r*(conc[0,k]+conc[0,k+1]+conc[2,k]-2*conc[1,k])
    
        for i in np.arange(2,xN-2):
            b[i-1] = 2*conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
            
            b[-1] = 2*conc[xN-2,k]+r*(conc[xN-1,k+1]+conc[xN-3,k]-2*conc[xN-2,k]+conc[xN-1,k])
    
        c = np.linalg.solve(A,b)
     
        for i in np.arange(0,xN-2):
            conc[i+1,k+1] = c[i]


    """ Plot solutions for all time """
    for i in np.arange(0,len(t),1):
       
        plt.plot(x,conc[:,i],'b',linewidth=3)
        plt.plot(-x,conc[:,i],'b',linewidth=3)
        plt.ylim((conc[-1,1],conc[0,1]))
        plt.pause(0.01)
        plt.clf()

    return conc


def crankydiffusionNOPLOT(x_start,x_stop,t_stop,bndryA,bndryB,iniconds,D,xN,tN):
        
    """ Set-up """
    # Step sizes
    h = (x_stop-x_start)/(xN-1)
    dt = (t_stop)/(tN-1) 
    r = (D*dt)/(h**2)
    
    # space and time vecotors, used to fill bndrys and iniconds
    x = np.linspace(x_start,x_stop,xN)
    
    # initialise solution matrix with given conditions
    conc = np.zeros((xN,tN))
    conc[0,:] = bndryA(x)
    conc[-1,:] = bndryB(x)
    conc[:,0] = iniconds(x)

    # initialise 'A' matrix (unknowns) to be solved at each time step 
    A = np.zeros((xN-2,xN-2))
    
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
    
    # set up b column vector (knowns)
    b = np.zeros((xN-2,1))
    
    
    """ Fill conc with approximations """
    for k in np.arange(0,tN-1):
    
        b[0] = 2*conc[1,k]+r*(conc[0,k]+conc[0,k+1]+conc[2,k]-2*conc[1,k])
    
        for i in np.arange(2,xN-2):
            b[i-1] = 2*conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
            
            b[-1] = 2*conc[xN-2,k]+r*(conc[xN-1,k+1]+conc[xN-3,k]-2*conc[xN-2,k]+conc[xN-1,k])
    
        c = np.linalg.solve(A,b)
     
        for i in np.arange(0,xN-2):
            conc[i+1,k+1] = c[i]


    """ Plot solutions for all time """
    

    return conc


