#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Section 19.4.2                                                     *
%*     filename: ch19pr04.m                                               *
%*     program listing number: 19.4                                       *
%*                                                                        *
%*     This program solves the Thomson problem  using genetic algorithm.  *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Improved by Alex Skinner.                                          *
%*     Last modification:  04/02/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# system parameters
N=5      # number of charges
R=1.0    # radius of the sphere

def U(theta,phi,r):  # potential energy (fitness function)
    u=0.0
    X=np.sin(theta)*np.cos(phi)
    Y=np.sin(theta)*np.sin(phi)    
    Z=np.cos(theta)    
    for i in range(0,N):
        for j in range(i+1,N):
            r12=np.sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2+(Z[i]-Z[j])**2)
            u=u+1.0/(r12*r)
    return u
    
# parameters for genetic algorithm
NP=2**12          # population size
MP=np.int(NP/2)   # half of the population
max_count=200;    # stopping condition
max_generation=10000   # maximum generation

# parameter space (genes)
th_max=np.pi; th_min=0.0
ph_max=2.0*np.pi; ph_min=0.0

# initial parameters (genes)
th=np.random.rand(N,NP)*(th_max-th_min)+th_min
ph=np.random.rand(N,NP)*(ph_max-ph_min)+ph_min
th[0,:]=0.0
ph[0,:]=0.0
ph[1,:]=0.0
th[2,0]=th_max
ph[2,1]=ph_max

# allocate arrays
ip=np.zeros(MP)
F=np.zeros(NP)
ths=np.zeros((N,NP))
phs=np.zeros((N,NP))

# reset the counters
found=False
count=0
generation=0
Fmin=np.finfo(np.float64()).max  # some large number 

# genetic evolution begins here
while not(found):
    
    generation=generation+1
    if generation > max_generation:  # too many generations
        print('Max generation reached.  Terminated.')
        break
    
    # eveluate the initial fitness
    for i in range(0,NP):
        F[i]=U(th[:,i],ph[:,i],R)
 
    # sort the population based on thier fitness
    IX=np.argsort(F)
    ths[:,:]=th[:,IX]; phs[:,:]=ph[:,IX] 
     
    # check if converged
    if F[IX[0]]>=Fmin:
        count=count+1
        if count > max_count:  # exceed the waiting time limit
            found=True
            break
    else:
        count=0   # reset wating time counter
        Fmin=F[IX[0]]
        print('New low found: generation={0:d}, fitness={1:.15f}'.format(generation,Fmin))
    
    # find mating pairs from the survived population
    ip=np.random.permutation(MP)
    
    # generate offsprings       
    for k in range(0,MP,2):
        # g ranges from -1 to 1
        g=-2.0*np.random.rand(2)+1.0
        i=ip[k]
        j=ip[k+1]
        ths[:,k+MP]  = ths[:,i]+g[0]*(ths[:,j]-ths[:,i])
        ths[:,k+1+MP]= ths[:,j]+g[1]*(ths[:,i]-ths[:,j])
        phs[:,k+MP]  = phs[:,i]+g[0]*(phs[:,j]-phs[:,i])
        phs[:,k+1+MP]= phs[:,j]+g[1]*(phs[:,i]-phs[:,j])

    # mutation rate              
    if np.mod(generation,10)==0:
        mutation=0.8*np.random.rand()
    else:
        mutation=0.3*np.random.rand()
        
    # mutation  
    rm=np.random.rand(NP)    
    for i in range(2,NP):
        if rm[i]<mutation:
            ths[:,i]=np.random.rand(N)*(th_max-th_min)+th_min
            phs[:,i]=np.random.rand(N)*(ph_max-ph_min)+ph_min
            ths[0,i]=0.0
            phs[0:1,i]=0.0
    
    # store new population
    th[:,:]=ths[:,:]; ph[:,:]=phs[:,:]

print('Lowest Energy={0:.15f}'.format(Fmin))
for i in range(0,N):
    print('theta={0:f}, phi={1:f}'.format(ths[i,0],phs[i,0]))

X=np.zeros(N)
Y=np.zeros(N)
Z=np.zeros(N)

for i in range(0,N):
    X[i]=R*np.sin(ths[i,0])*np.cos(phs[i,0])
    Y[i]=R*np.sin(ths[i,0])*np.sin(phs[i,0])
    Z[i]=R*np.cos(ths[i,0])

# plot 3D phase trajectory
plt.close('all')
fig=plt.figure()
ax=fig.gca(projection='3d')
plt.plot(X, Y, Z,'ob')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.paspect=[1.0,1.0,1.0]
ax.set_xlim3d([-1.0,1.0])
ax.set_ylim3d([-1.0,1.0])
ax.set_zlim3d([-1.0,1.0])

bond=10.0
for i in range(0,N):
    for j in range(i+1,N):
        d=np.sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2+(Z[i]-Z[j])**2)
        bond = np.min([bond,d])

bond=bond*1.25
for i in range(0,N):
    for j in range(i+1,N):
        d=np.sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2+(Z[i]-Z[j])**2)
        if d <= bond:
            ax.plot([X[i],X[j]],[Y[i],Y[j]],[Z[i],Z[j]],'-b')
            plt.pause(0.0001)
plt.show()