#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 16:14:16 2020

@author: WillMcD
"""
import matplotlib.pyplot as plt
import numpy as np
import os #,sys


class Spark:
    def __init__(self,location,fire_size):
        self.location = location
        self.fired = False
        self.time = None
        self.fire_size = fire_size


""" Set-up """
# Step sizes
x_start = -10 # bound for the LHS of the x axis
x_stop = 10 # bound for the RHS of the x axis
t_stop = 2 # How long you want to the simulation to run for

D = 1
a = 0.005


h = 0.01
dt = 0.01



# space and time vecotors, used to fill bndrys and iniconds
x = np.arange(x_start,x_stop,h)
t = np.arange(0,t_stop,dt)

spark_steps = 100


def pde_implicit(D,x,t,x_start,x_stop,t_stop,spark_steps):
    
    xN = ((x_stop-x_start)/h) + 1
    
    r = (D*dt)/((h**2))
    
    # =============================================================================
    #  setting up lists/counters for later on
    calcium_conc_list = []
    spark_list = []
    diff_list = []
    
    # setting up each spark
    for i in range((x_stop-x_start + 1)*2):
        spark_list.append(Spark(i,1/a))
    # =============================================================================
    
    
    conc = np.zeros((len(x),len(t)))
    conc[0,:] = 0
    
    # initialise 'A' matrix (unknowns) to be solved at each time step 
    A = np.zeros((len(x)-2,len(x)-2))
    
    conc[-1,:] = 0
    conc[int(len(x)/2)+1,0] = 600
    
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
            if i >= spark_steps and i%spark_steps == 0:
                location = int(i/spark_steps)
                spark = spark_list[location]
                if spark.fired is False and c[i] >= 1:
                    spark.location = location
                    print('Firing at location = ' + str(spark.location-int((x_stop-x_start)/2)) + '    at time ' + str(k))
                    c[i] += spark.fire_size
                    spark.fired = True
                    spark.time = k
            conc[i+1,k+1] = c[i]
        
        calcium_conc_list.append(total_conc)
        conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]
        conc[-1,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]
        
     
    for spark in spark_list:
        half_point = int((x_stop-x_start)/2)
        RHS_start = half_point + 2
        if spark.location >= RHS_start and spark.time is not None:
            difference = spark.time - spark_list[spark.location-1].time
            diff_list.append(difference)
    
    return [diff_list, conc]
    
    
def plot(conc, x, t):
    
    plot_label_count_numeric = 0
    
    fig, ax = plt.subplots()
    filename = 'foo0000.jpg'
    plt.savefig(filename,format='jpg')
    c = 0
    
    for i in np.arange(0,len(t),1):
        ax.plot(x,conc[:,i],'b',linewidth=3, label = str(plot_label_count_numeric))
        ax.set_ylim(0,8)
        
        print(c)
        filename = 'foo' + str(c+1).zfill(4) + '.jpg'
        plt.savefig(filename, format='jpg')

        plt.pause(0.01)
        line = [line for line in ax.lines if line.get_label() == str(plot_label_count_numeric)][0]
        ax.lines.remove(line)
        plot_label_count_numeric += 1
        c = c + 1
    
    os.system("ffmpeg -y -i foo%04d.jpg implicit.m4v")
    os.system("rm -f *.jpg")
    
    return

result = []

for D in [0.5,1,1.5,2,2.5,3]:
    matrix = pde_implicit(D,x,t,x_start,x_stop,t_stop,spark_steps)
    result.append(matrix)

#plot(result[1][1],x,t)

time_result_list = []

for i in result:
    diff_list = i[0]
    difference = None
    for k in diff_list:
        k_1 = None
        k_2 = None
        k_3 = None
        if k_3 is not None:
            if k_2 is not None:
                if k_1 is not None:
                    if k == k_1 or k == k_2 or k == k_3:
                        difference = k
                        break
                    else:
                        k_3 = k_2
                        k_2 = k_1
                        k_1 = k
                else:
                    k_3 = k_2
                    k_2 = k_1
                    k_1 = k
            else:
                k_2 = k_1
                k_1 = k
        else:
            k_1 = k
    if difference is None:
        difference = k_3
    time_result_list.append(difference)
        


    