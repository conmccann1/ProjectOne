import numpy as np
import math
import pylab as pl
import time

start = time.time()

n = 1000
alpha = 1/np.sqrt(2*np.pi*np.e)

x = np.linspace(-3,3,100)
c = np.linspace(0,4,100)
t = np.linspace(0.1, 10, n)


def c_func(x,t):
    return np.exp((-x**2/(4*t)))/math.sqrt(4*np.pi*(alpha**2)*t)


fig, ax = pl.subplots()

ax.plot(x,np.ones(100),linestyle='dashed', color = 'red', linewidth = 0.5)
ax.plot(np.ones(100),c,linestyle='dotted', color = 'green', linewidth = 0.5)
ax.plot(np.negative(np.ones(100)),c,linestyle='dotted', color = 'green', linewidth = 0.5)
ax.set_ylabel(r'$\bar{c}$', fontsize = 20)
ax.set_xlabel(r'$\bar{x}$', fontsize = 20)

c = 0

for i in t:
      
    curve = []

    for j in x:
        curve.append(c_func(j,i))
    
    if c == 0 or c == 50:
        ax.plot(x,curve, color = 'black', linewidth = 0.5, label = str(c))
        pl.pause(0.01)

    elif c > 150:
        break
    else:
        print('c = ' + str(c))
        ax.plot(x,curve, color = 'black', linewidth = 0.5, label = str(c))
        pl.pause(0.01)
        line = [line for line in ax.lines if line.get_label()==str(c)][0]
        ax.lines.remove(line)
    c += 1


print('end: ' + str(time.time()-start))