# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:22:48 2020

@author: Conor
"""

#from CrankNicolsonNeumann import crank_nicolson_neu
#import CrankNicolsonDirichlet as cnd
import matplotlib.pyplot as plt
import numpy as np
import os #,sys
import UsefulFunctions as uf
from waterfall import waterfall


class Spark:
    def __init__(self,location,fire_size):
        self.location = location
        self.fired = False
        self.time = None
        self.fire_size = fire_size



""" Set-up """
# Step sizes
x_start = -5 # bound for the LHS of the x axis
x_stop = 5 # bound for the RHS of the x axis
#xN = 300
#tN = 300
t_stop = 1.8 # How long you want to the simulation to run for
#dt = (t_stop)/(tN-1)
D = 1
a = 0.005


h = 0.001
dt = 0.001

xN = ((x_stop-x_start)/h) + 1
tN = ((t_stop)/dt) + 1

#h = (x_stop-x_start)/(xN-1)
r = (D*dt)/((h**2))

# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop,h)
t = np.arange(0,t_stop,dt)


# =============================================================================
#  setting up lists/counters for later on
calcium_conc_list = []
spark_list = []
diff_list = []
progress = 0

# setting up each spark
for i in range((x_stop-x_start + 1)*2):
    spark_list.append(Spark(i,1/a))
# =============================================================================


conc = np.zeros((len(x),len(t)))
#    conc[-1,:] = 2
#    conc[:,0] = 1
conc[0,:] = 0

# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((len(x)-2,len(x)-2))

conc[-1,:] = 0
conc[int(len(x)/2)+1,0] = 6000
#A[0,0] = r*(2/3)+2
#A[0,1] = -r*(2/3)
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
  
for k in np.arange(0,len(t)-1,1):
    
    total_conc = 0
    
    for i in np.arange(1,len(x)-1):
        b[i-1] = conc[i,k]
        
        
    c = np.linalg.solve(A,b)

   
    
    for i in np.arange(0,len(x)-2):
        total_conc += c[i]
        if i >= 1000 and i%1000 == 0:
            location = int(i/1000)
            spark = spark_list[location]
            if spark.fired is False and c[i] >= 1:
                spark.location = location
                print('Firing at location = ' + str(spark.location-int((x_stop-x_start)/2)) + '    at time ' + str(k))
                c[i] += spark.fire_size
                spark.fired = True
                spark.time = k
        conc[i+1,k+1] = c[i]
    
    calcium_conc_list.append(total_conc)
    # print('Total concentration = ' + str(total_conc))
    conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]
    conc[-1,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]
    
# =============================================================================
#     percent = len(np.arange(0,len(t)-1,1))/20
#     
#     if k%percent == 0:
#         if progress == 0:
#             print('Calculation progress: 20%')
#             progress += 1
#         if progress == 1:
#             print('Calculation progress: 40%')
#             progress += 1
#         if progress == 2:
#             print('Calculation progress: 60%')
#             progress += 1
#         if progress == 3:
#             print('Calculation progress: 80%')
#             progress += 1
#         if progress == 4:
#             print('Calculation progress: 100%')
# =============================================================================


for spark in spark_list:
    half_point = int((x_stop-x_start)/2)
    RHS_start = half_point + 2
    if spark.location >= RHS_start and spark.time is not None:
        # print('spark location: ' + str(spark.location))
        difference = spark.time - spark_list[spark.location-1].time
        diff_list.append(difference)


err = np.zeros((len(x),len(t)))
errSumX = np.zeros((len(x),1))


plot_label_count_numeric = 0
# plot_label_count_analytic = 0

fig, ax = plt.subplots()

plt.plot(np.ones(100),np.arange(0,10,0.1),'k--', linewidth = 0.5)
plt.plot(-1*np.ones(100),np.arange(0,10,0.1),'k--', linewidth = 0.5)
plt.plot(x,np.ones(len(x)),'k--', linewidth = 0.5)
plt.xlim(-5,5)
plt.ylim(0,8)
plt.xlabel(r'$\bar{x}$', fontsize = 20)
h = plt.ylabel(r'$\bar{c}$', fontsize = 20)
h.set_rotation(0)


#filename = 'foo0000.jpg'
#plt.savefig(filename,format='jpg')
c = 0
for i in np.arange(0,len(t),1):
    err[:,i] = abs(conc[:,i]-uf.cBar(x,t[i])).T
    ax.plot(x,conc[:,i],'b',linewidth=3, label = str(plot_label_count_numeric))
    # ax.plot(x[::20],uf.cBar(x[::20],t[i]),'rx',mew=5,ms=5, label = str(plot_label_count_analytic))
    ax.set_ylim(0,8)
    
    if i == spark_list[2].time-3:
        plt.savefig('beforespark.jpg', format='jpg')
        print('time before = ' +str(i))
    if i == spark_list[2].time+3:
        plt.savefig('afterspark.jpg', format='jpg')
        print('time after = ' +str(i))
#    print(c)
#    filename = 'foo' + str(c+1).zfill(4) + '.jpg'
#    plt.savefig(filename, format='jpg')
    
##    plt.ylim((0,cmax))
#    plt.xlim((x_start,3))
    plt.pause(0.001)
    line = [line for line in ax.lines if line.get_label() == str(plot_label_count_numeric)][0]
    ax.lines.remove(line)
    # line = [line for line in ax.lines if line.get_label() == str(plot_label_count_analytic)][0]
    # ax.lines.remove(line)
    plot_label_count_numeric += 1
    # plot_label_count_analytic += 1
#    c = c + 1


#os.system("ffmpeg -y -i foo%04d.jpg implicit.m4v")
#os.system("rm -f *.jpg")

#waterfall(conc,t)

# =============================================================================
# errSumX = np.sum(err[:,1:],axis=1)
# errSumT = np.sum(err[:,1:],axis=0)
# 
# plt.figure()
# plt.subplot(211)
# plt.title('Error for CN+I')
# plt.plot(x,errSumX,linewidth=2)
# plt.xlabel('x')
# plt.ylabel('sum of errors')
# plt.subplot(212)
# plt.plot(t[1:],errSumT,linewidth=2)
# plt.xlabel('t')
# =============================================================================
