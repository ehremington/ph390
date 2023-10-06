#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%**************************************************************************
%*     Example 19.2                                                       *
%*     filename: ch19pr01.py                                              *
%*     program listing number: 19.1                                       *
%*                                                                        *
%*     This program attempt to find a grobal minimum of a fitness         *
%*     function U(x) using simulated annealing.                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/01/2017.                                    *
%**************************************************************************
"""
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# system setting
L=241    # size of the descretized configuration space
dx=0.1   # grid size
xmin=0.0 # lower bound of the configuration space
xmax=(L-1)*dx  # upper boound
x=np.linspace(xmin,xmax,L)

def U(x):  # fitness function
    return -3.*np.cos(2.*(x-xmax/2.))+0.2*(x-xmax/2.)**2+3.0

# graphics parameter
movie=False   # show animation

NP=2**3   # population size
p=np.random.rand(NP)*(xmax-xmin)+xmin # initial configuration
T=5.0 # initial temperature
E=U(p)  # initial fitness
dp=0.1  # Metropolis max step length

NT=200   # Total cooling steps
NS=2000  # Total thermalization steps
RT=0.98  # Cooling rate

# allocate arrayss
temp=np.zeros(NT)
Emin=np.zeros(NT)

if movie:
    plt.figure(figsize=(6,5))
    
for k in range(0,NT): # loop over time step
    temp[k]=T

    for i in range(0,NS): # Thermalization Loop (Metropolis)
        
        for j in range(0,NP):  # Loop over population

            found=False
            while not(found):
                p0=p[j]+np.random.choice([-1,1])*dp
                E0=U(p0)
                if np.exp(-(E0-E[j])/T) > np.random.rand(1):
                    found=True
                    if p0<xmin or p0>xmax:
                        p[j]=np.random.rand(1)*(xmax-xmin)+xmin
                        E[j]=U(p[j])
                    else:
                        p[j]=p0
                        E[j]=E0

        if movie: # update movie
            plt.clf()
            plt.plot(x,U(x))
            plt.xlim([xmin,xmax])
            plt.ylim([-1,19])

            for j in range(0,NP):
                circle=plt.Circle((p[j],U(p[j])),0.25,fc='b')
                plt.gca().add_patch(circle)
                plt.pause(0.0001)

    Emin[k]=np.min(E)
    print('T={0:f},   Emin={1:f}'.format(T,Emin[k]))
    
    T=T*0.98  # Exponential cooling schedule

print('Final Temperature = {0:f}'.format(temp[k]))
print('Final Lowest Fitness = {0:f}'.format(Emin[k]))
    
if not(movie):
    plt.figure(figsize=(6,5))
    plt.plot(x,U(x))
    plt.xlim([xmin,xmax])
    plt.ylim([-1,19])

    for j in range(0,NP):
        circle=plt.Circle((p[j],U(p[j])),0.25,fc='b')
        plt.gca().add_patch(circle)
        plt.pause(0.0001)

plt.xlabel('x',fontsize=14)
plt.ylabel('Fitness',fontsize=14)
plt.show()
    
plt.figure(figsize=(6,5))
t=np.linspace(1,NT,NT)
plt.plot(t,temp,'-k',label='temperature')
plt.plot(t,Emin,'-r',label='energy')
plt.xlabel('Steps',fontsize=14)
plt.ylabel('Energy',fontsize=14)
plt.legend(loc=1)
plt.show()
    
 
