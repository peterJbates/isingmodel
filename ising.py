import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from math import exp
from time import time

#----------------
#Global Variables
#----------------

start_time = time()
n = 20                 #used to create n x n domain
N = n**2               #number of atoms in domain
iterations = 10**5

#creates a normalized distribution of temperature points around the mean_temp
mean_temp = 2.27
temps = np.random.normal(mean_temp, .64, N)
temps = temps[(temps>1) & (temps<4)]
n1 = 1/(iterations*N)
nt = np.size(temps)
magnetization = np.zeros(nt)
energy = np.zeros(nt)

specific_temps = [1,1.5,2,2.27,2.5,3] #temperatures used to create domain images


#-------------------------------
#Functions Utilized by Main Loop
#-------------------------------

'''Creates a random array of dipole states'''
def initialize(n_sites):
    state = np.zeros((n_sites, n_sites))

    for i in range(n_sites):
        for j in range(n_sites):
            if np.random.random() < 0.5:
                state[i][j] = 1
            else:
                state[i][j] = -1
    return state

'''Computes dU of a flipping a dipole.
This is a Monte Carlo algorithm using
the Metropolis Algorithm. It has peridic
boundary conditions.'''
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

'''Computes overall magnetism of a state'''
def Magnetization(state):
    mag = np.sum(state)
    return mag

'''Computes overall energy of a state'''
def Energy(state):
    global n
    energy = 0
    for i in range(len(state)):
        for j in range(len(state)):
            S = state[i,j]
            nb = state[(i+1)%n, j] + state[i, (j+1)%n] + state[(i-1)%n, j] + state[i,(j-1)%n]
            energy += -nb*S
    return energy/4


#---------
#Main Loop
#---------

for m in range(len(temps)):
    E1 = M1 = 0
    T = temps[m]
    state = initialize(n)  #calls the initialize function

    for x in range(iterations):
        i = int(np.random.random()*n)
        j = int(np.random.random()*n)
        Ediff = deltaU(i,j,n,state) #calls the deltaU function for random dipole
        if Ediff <= 0:
            state[i,j] *= -1
        else:
            if np.random.random() < exp(-Ediff/T):
                state[i,j] *= -1
        Ene = Energy(state)
        Mag = Magnetization(state)

        E1 += Ene
        M1 += Mag

    magnetization[m] = n1*M1
    energy[m] = n1*E1


print("---%s seconds---" %(time()-start_time))  #prints main loop run time


#---------------------
#Graphs and Animations
#---------------------

'''Creates a list of domain equilibrium states at specific temperatures.
This list will be used to create images.'''
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))
for m in range(len(specific_temps)):
    T = specific_temps[m]
    state = initialize(n)

    for x in range(iterations):
        i = int(np.random.random()*n)
        j = int(np.random.random()*n)
        Ediff = deltaU(i,j,n,state)
        if Ediff <= 0:
            state[i,j] *= -1
        else:
            if np.random.random() < exp(-Ediff/T):
                state[i,j] *= -1
    if m < 3:
        axs[0,m].set_title('T=%s' %(specific_temps[m]))
        axs[0,m].imshow(state, cmap='tab20c')
        axs[0,m].set_xticks([])
        axs[0,m].set_yticks([])
    else:
        axs[1,m-3].set_title('T=%s' %(specific_temps[m]))
        axs[1,m-3].imshow(state, cmap='tab20c')
        axs[1,m-3].set_xticks([])
        axs[1,m-3].set_yticks([])

plt.savefig('eqstates.jpg', format='jpg')
plt.show()
plt.close(fig)


'''This Uses the results from the main loop to create
graphs for magnetization/spin and energy/spin.'''
fig = plt.figure(figsize = (15,5))

sub = fig.add_subplot(1,2,1)
plt.plot(temps, abs(magnetization), 'o', markerfacecolor= "blue", markersize= 3)
plt.xlabel("Temperature (T)", fontsize=12)
plt.ylabel("Magnetization per site", fontsize=12)

sub = fig.add_subplot(1,2,2)
plt.plot(temps, energy, 'o', markerfacecolor = 'green', markersize = 3)
plt.xlabel("Temperature (T)", fontsize=12)
plt.ylabel("Energy per site", fontsize=12)

plt.savefig("magnetizationAndEnergy.jpg", format='jpg')
plt.show()
plt.close(fig)
