import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from math import exp
from time import time

#------------------
#Functions Utilized
#------------------

#Creates a random array of dipole states
def initialize(n_sites):
    state = np.zeros((n_sites, n_sites), dtype=np.int8)

    for i in range(n_sites):
        for j in range(n_sites):
            if np.random.random() < 0.5:
                state[i][j] = 1
            else:
                state[i][j] = -1
    return state


#Computes dU of flipping a dipole
def deltaU(i,j,n_sites,state):
    m = n_sites-1    #max row/column entry

    if i == 0:
        top = state[m,j]
    else:
        top = state[i-1,j]

    if i == m:
        bottom = state[0,j]
    else:
        bottom = state[i+1,j]

    if j == 0:
        left = state[i,m]
    else:
        left = state[i, j-1]

    if j == m:
        right = state[i,0]
    else:
        right = state[i,j+1]

    dU = 2*state[i,j]*(top+bottom+left+right)
    return dU

def Magnetization(state):
    mag = np.sum(state)
    return mag


#---------
#Main Loop
#---------

start_time = time()


n = 10
iterations = 100 * n**2     #change back to 100 * n**2
mean_temp = 2.27
temps = np.random.normal(mean_temp, .64, n**2)
temps = temps[(temps>1) & (temps<4)]
nt = np.size(temps)
mag = np.zeros(nt)


for m in range(len(temps)):
    T = temps[m]
    state = initialize(n)

    for x in range(iterations*100):
        i = int(np.random.random()*n)
        j = int(np.random.random()*n)
        Ediff = deltaU(i,j,n,state)
        if Ediff <= 0:
            state[i,j] *= -1
        else:
            if np.random.random() < exp(-Ediff/T):
                state[i,j] *= -1

    mag[m] = Magnetization(state)/n**2

print(temps, mag)
print("---%s seconds---" %(time()-start_time))


#---------------------
#Graphs and Animations
#---------------------

fig = plt.figure()

sub = fig.add_subplot(2,2,2)
plt.plot(temps, abs(mag), 'o', color='blue')
plt.xlabel("Temperature (T)", fontsize=12)
plt.ylabel("Magnetization", fontsize=12)
plt.show()
