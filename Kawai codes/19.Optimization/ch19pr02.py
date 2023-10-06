#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 19.3                                                       *
%*     filename: ch19pr02.py                                              *
%*     program listing number: 19.2                                       *
%*                                                                        *
%*     This program findw a grobal minimum of a fitness function U(x)     *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/01/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

#system configuration
L=121          # size of the descretized configuration space
dx=0.1         # grid size
xmin=0.0       # lower bound of the configuration space
xmax=(L-1)*dx  # upper bound
x=np.linspace(xmin,xmax,L)

def U(x): # fitness function
    return np.cos(5.0*x)-2.0*np.sin(3.5*x)+0.5*np.cos(x+0.5)+4.0 

# parameter sof genetic algorithm
N=2**8            # size of population
M=np.int(N/2)
p=np.random.rand(N)*xmax  # initial population

max_count=100     # waiting time
mutaion=0.1       # mutation rate (fixed)
max_generation=10000
bestfit=np.zeros(max_generation)

# initial fitness
f=U(p)  # fitness
ix=np.argsort(f)   # sorting population basd on their fitness
fmin=f[ix[0]] 
q=p[ix]
bestfit[0]=fmin

# plot initial population
plt.figure(figsize=(6,6))
plt.plot(x,U(x))
plt.xlim([xmin, xmax])
plt.ylim([-1, 11])
plt.axis('equal')

for i in range(0,N):
    circle=plt.Circle((p[i],U(p[i])),0.2,facecolor='none', edgecolor='b')
    plt.gca().add_patch(circle)
plt.pause(0.0001)

count=0;
found=False
ip=np.zeros(M)
generation=0;

while not(found):
    generation=generation+1
    if generation > max_generation:
        print('Max generation reached.  Terminated.')
        break
    
    ip=np.random.permutation(M)  
    
    # generate offspring
    for k in range(0,M,2):

        i=ip[k]
        j=ip[k+1]
        g=np.random.rand(2)
        # replace the deads with the new borns
        q[k+M]= q[i]+g[0]*(q[j]-q[i]+0.01)     # inheritance
        q[k+1+M]= q[j]+g[1]*(q[i]-q[j]+0.01)  # inheritance
    
    # mutation
    if np.mod(generation,10)==0:
        mutation=0.8
    else:
        mutation=0.1
    
    if np.random.rand(1)<mutaion:
        i=np.random.randint(0,N)
        q[i]=np.random.rand(1)*(xmax-xmin)

    p=q  # store new population
    
    # evaluate fitness
    f=U(p)
    plt.clf()
    plt.plot(x,U(x))
    plt.xlim([xmin, xmax])
    plt.ylim([-1, 11])

    for i in range(0,N):
        circle=plt.Circle((p[i],U(p[i])),0.2,facecolor='none', edgecolor='b')
        plt.gca().add_patch(circle)
    plt.pause(0.0001)

    #  sort the population
    ix=np.argsort(f)
    q=p[ix]
    
    #check if converged
    if f[ix[0]]>=fmin :
        count=count+1
        if count > max_count:
            found=True  # waited long enough

    else:
        count=0  # reset wating couter
        fmin=f[ix[0]]
        print('New low found: generation={0:d} fitness={1:.15f}'.format(generation,fmin))

    bestfit[generation]=fmin

print('optimal x={0:f},   U(x)={1:f}'.format(q[0],f[ix[0]]))

plt.figure(figsize=(6,5))
t=np.linspace(0,generation,generation+1)
plt.plot(t,bestfit[0:generation+1],'-or',mfc=None)
plt.xlabel('Generation',fontsize=14)
plt.ylabel('Best Finess',fontsize=14)
plt.show()
