import matplotlib.pyplot as plt
import numpy as np

x_start = 0
x_stop = 12

t_stop = 1
D = 1
xN = int(50*x_stop)
tN = int(200*t_stop)
#h = 0.04291845
#dt = 0.00505050505
a = 0.053 #(1/((2*np.pi*np.e)**0.5))

plot_label_count_numeric = 0
plot_label_count_analytic = 0

""" Set-up """
# setting up list for sparks

location_list = []
fired = []
for i in range(x_stop):
    fired.append(False)

# Step sizes
h = (x_stop-x_start)/(xN-1)
dt = (t_stop)/(tN-1)
r = (D*dt)/(h**2)

# space and time vecotors, used to fill bndrys and iniconds
x = np.linspace(x_start,x_stop,xN)
# print(x)
t = np.linspace(0,t_stop,tN)
#t = np.arange(0,t_stop,dt)

# initialise solution matrix with given conditions
conc = np.zeros((xN,tN))
conc[-1,:] = 0
conc[0,0] = 1/h


# initialise 'A' matrix (unknowns) to be solved at each time step 
A = np.zeros((xN-2,xN-2))

A[0,0] = (r*(2/3)+2)
A[0,1] = -r*(2/3)

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

    b[0] = 2*conc[1,k]+r*(conc[0,k]-2*conc[1,k]+conc[2,k])

    for i in np.arange(2,xN-2):
        b[i-1] = 2*conc[i,k]+r*(conc[i-1,k]-2*conc[i,k]+conc[i+1,k])
        
    b[-1] = 2*conc[xN-2,k]+r*(+conc[xN-3,k]-2*conc[xN-2,k])
    b[-1] = 2*conc[xN-2,k]+r*(conc[xN-1,k+1]+conc[xN-3,k]-2*conc[xN-2,k]+conc[xN-1,k])

    c = np.linalg.solve(A,b)
 
    for i in np.arange(0,xN-2):
# =============================================================================
#         if i != 0 and i%50 == 0:
#             location = int(i/50)
#             if fired[location] is False and c[i] >= 1:
#                 print('Firing at location = ' + str(location))
#                 c[i] += 1/h
#                 fired[location-1] = True
#                 fired[location] = True
# =============================================================================

        conc[i+1,k+1] = c[i]
    
    conc[0,k+1] = (4/3)*conc[1,k+1]-(1/3)*conc[2,k+1]


def cBar(x,t):
    return np.exp(-(x**2) / (4*t)) / ( ((4*t*np.pi*(a**2)))**0.5 )


fig, ax = plt.subplots()
ax.plot(np.ones(40),np.arange(0,4,0.1),'--')
ax.plot(-1*np.ones(40),np.arange(0,4,0.1),'--')
ax.plot(x,np.ones(len(x)),'--')
ax.plot(-x,np.ones(len(x)),'--')

""" Plot solutions for all time """
for i in np.arange(0,len(t),1):
    ax.plot(x[::5],cBar(x[::5],t[i]),'xb',mew=5, ms=5, label = 'R' + str(plot_label_count_analytic))
    ax.plot(-x[::5],cBar(x[::5],t[i]),'xb',mew=5, ms=5, label = 'L' + str(plot_label_count_analytic))
    ax.plot(x,conc[:,i],'r',linewidth=2, label = 'R' + str(plot_label_count_numeric))
    ax.plot(-x,conc[:,i],'r',linewidth=2, label = 'L' + str(plot_label_count_numeric))
    ax.set_ylim((0,7))
    ax.set_xlim((-12,12))
    plt.pause(0.01)
    line = [line for line in ax.lines if line.get_label()=='R' + str(plot_label_count_analytic)][0]
    ax.lines.remove(line)
    line = [line for line in ax.lines if line.get_label()=='R' + str(plot_label_count_numeric)][0]
    ax.lines.remove(line)
    line = [line for line in ax.lines if line.get_label()=='L' + str(plot_label_count_analytic)][0]
    ax.lines.remove(line)
    line = [line for line in ax.lines if line.get_label()=='L' + str(plot_label_count_numeric)][0]
    ax.lines.remove(line)
    plot_label_count_numeric += 1
    plot_label_count_analytic += 1

