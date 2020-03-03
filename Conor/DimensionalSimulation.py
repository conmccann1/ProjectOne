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
from matplotlib import cm
a = (1/((2*np.pi*np.e)**0.5))

""" Set-up """
# Step sizes
x_start = -10
x_stop = 10
#xN = 300
#tN = 300
t_stop = 3
#dt = (t_stop)/(tN-1)
D = 1
d = 2
#D = 30


h = 0.01
dt = 0.01

xN = ((x_stop-x_start)/h) + 1
tN = ((t_stop)/dt) + 1

#h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)
r = (D*dt)/((h**2))

# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop,h)
t = np.arange(0,t_stop,dt)
#print(t)

siteloc = np.arange(1/h,2*x_stop/h,1/h)


conc = np.zeros((len(x),len(t)))
#    conc[-1,:] = 2
#    conc[:,0] = 1
conc[0,:] = 0

# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((len(x)-2,len(x)-2))

conc[-1,:] = 0
conc[int(len(x)/2),0] = 430#430

A[0,0] = 1 + (2/3)*r
A[0,1] = -(2/3)*r

A[-1,-1] = 1 + (2/3)*r
A[-1,-2] = -(2/3)*r
for i in np.arange(1,len(x)-3,1):
    A[i,i] = 2*(r)+1
    A[i,i-1] = -r
    if i == xN-3:
        continue
    A[i,i+1] = -r
    
# set up b column vector (knowns)
b = np.zeros((len(x)-2,1))

a = 0.05
siteL = 1
siteR = 1
sparktimes = np.zeros(x_stop*2)

sparktimes = np.zeros(0)

for k in np.arange(0,len(t)-1,1):
    
        
        
    for i in np.arange(int(1/h),int((2*x_stop/h)/2),int(1/(h))):
#        print(i)
        if abs(abs(x[i])-siteL) <= 0.0000001 and conc[i,k] >= 1.0:
#            print(conc[i,k])
            conc[i,k]+=430
            siteL += 1
            print('SPARK')
            print('t= %f' %t[k])
            print('t= %f' %k)
            print('site = %f' %x[i])
            print('site = %f' %i)
            sparktimes = np.append([sparktimes],[t[k]])
            print(sparktimes)
#            print('\n')
#            print(len(conc[0:i,k]))
        
    for i in np.arange(int((2*x_stop/h)/2),int((2*x_stop/h)),int(1/(h))):
        
        if abs(x[i]-siteR) <= 0.0000001 and conc[i,k] >= 1.0:

#            print(conc[i,k])
#            print('SPARK')
            conc[i,k]+=430
            siteR += 1
            print('t= %f' %t[k])
            print('site = %f' %x[i])
            print('\n')
#            print(len(conc[0:i,k]))
            
        else:
            continue
        
    for i in np.arange(1,len(x)-1):
        b[i-1] = conc[i,k]
    c = np.linalg.solve(A,b)
    
    for i in np.arange(0,len(x)-2):
        conc[i+1,k+1] = c[i] 
    
    conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]
    conc[-1,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]

std = np.diff(sparktimes)
print(sparktimes)
print(std)
#sm1 = conc[400,:]
#s1 = conc[600,:]
#
#sdif = abs(sm1-s1)

#err = np.zeros((len(x),len(t)))
#errSumX = np.zeros((len(x),1))
D=30
d=2
x = x*2

t = t*((d**2)/D)
#std = std*((d**2)/D)
sparktimes = sparktimes*((d**2)/D)
std = np.diff(sparktimes)
conc = 0.14*conc
plt.figure()
#for i in np.arange(0,len(t),1):
##    err[:,i] = abs(conc[:,i]-uf.cBar(x,t[i])).T
#    plt.plot(x,conc[:,i],'b',linewidth=3)
##    plt.plot(x[::5],uf.cBar(x[::5],t[i]),'r',mew=5,ms=5)
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
##    plt.xlim((-5,5))
#    plt.pause(0.00000001)
#    plt.clf()
#[CONC,T] = np.meshgrid(conc,t)
#[x,t] = np.meshgrid(x,t)
levels = np.arange(conc.min(),conc.max(),0.1)
norm = cm.colors.Normalize(vmax=0.1*conc.max(),vmin=conc.min())
plt.contourf(x,t,conc.T,norm=norm,levels=levels)   
plt.ylim(max(t), min(t))
#errSumX = np.sum(err[:,1:],axis=1)
#errSumT = np.sum(err[:,1:],axis=0)
#
#plt.figure()
#plt.subplot(211)
#plt.title('Error for CN+I')
#plt.plot(x,errSumX,linewidth=2)
#plt.xlabel('x')
#plt.ylabel('sum of errors')
#plt.subplot(212)
#plt.plot(t[1:],errSumT,linewidth=2)
#plt.xlabel('t')
#plt.close('all')